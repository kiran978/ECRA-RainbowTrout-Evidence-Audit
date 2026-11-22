#!/usr/bin/env python3
"""Research-only external rainbow-trout QSAR benchmarking.

NOT a controlled release of ECRA primary-study records. The public secondary
source does not contain per-measurement primary-study reference IDs or human approvals.
"""
from pathlib import Path
import argparse, hashlib, json, math, re, warnings
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, Lipinski
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
import rdkit, sklearn
from sklearn.model_selection import GroupShuffleSplit, GroupKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import Ridge
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
RDLogger.DisableLog('rdApp.*')
warnings.filterwarnings('ignore', category=FutureWarning)

ROOT=Path(__file__).resolve().parent
PINPUT=ROOT/'input'/'public_RT_secondary_experimental.csv'
FROZEN=ROOT/'input'/'ECRA_frozen_169_HOLD.csv'
OUT=ROOT/'output';OUT.mkdir(exist_ok=True)
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
MOD=ROOT/'models';MOD.mkdir(exist_ok=True)

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def cas_valid(x):
 if not isinstance(x,str):return False
 s=x.strip()
 if not re.fullmatch(r'[0-9]{2,7}-[0-9]{2}-[0-9]',s):return False
 a,b,c=s.split('-');digits=a+b
 return sum(int(t)*(i+1) for i,t in enumerate(digits[::-1]))%10==int(c)

def metric(y,yp):
 return {'n':len(y),'MAE':float(mean_absolute_error(y,yp)), 'RMSE':float(np.sqrt(mean_squared_error(y,yp))), 'R2':float(r2_score(y,yp)) if len(y)>1 else None}

def descriptor(m):
 return [float(Descriptors.MolWt(m)),float(Crippen.MolLogP(m)),float(rdMolDescriptors.CalcTPSA(m)),
 float(Lipinski.NumHDonors(m)),float(Lipinski.NumHAcceptors(m)),float(Lipinski.NumRotatableBonds(m)),
 float(rdMolDescriptors.CalcNumRings(m)),float(rdMolDescriptors.CalcFractionCSP3(m)),
 float(m.GetNumHeavyAtoms()),float(rdMolDescriptors.CalcNumAromaticRings(m)),
 float(rdMolDescriptors.CalcNumHeteroatoms(m)),float(Chem.GetFormalCharge(m))]

def split_features(d):
 fg=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=512)
 desc=[];bits=[];fps=[]
 for sm in d['Canonical_Isomeric_SMILES']:
  m=Chem.MolFromSmiles(sm); desc.append(descriptor(m))
  fp=fg.GetFingerprint(m)
  v=np.zeros((512,),dtype=np.int8);DataStructs.ConvertToNumpyArray(fp,v)
  bits.append(v);fps.append(fp)
 return np.hstack([np.asarray(bits,dtype=np.float32),np.asarray(desc,dtype=np.float32)]),np.asarray(desc,dtype=np.float32),fps

def get_models(seed=20261122):
 return {
 'Ridge':make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),Ridge(alpha=75.0)),
 'RandomForest':RandomForestRegressor(n_estimators=180,max_features=0.45,min_samples_leaf=2,n_jobs=-1,random_state=seed),
 'ExtraTrees':ExtraTreesRegressor(n_estimators=220,max_features=0.55,min_samples_leaf=2,n_jobs=-1,random_state=seed)}

def scaffold_of(m,canonical):
 s=MurckoScaffoldSmiles(mol=m,includeChirality=False)
 return 'SCAFFOLD_'+s if s else 'ACYCLIC_'+canonical

def main():
 source=pd.read_csv(PINPUT,dtype={'CAS No.':str,'uid':str})
 frozen=pd.read_csv(FROZEN,dtype={'CAS':str,'Result_ID':str})
 frozen_cas=set(frozen.CAS.dropna().astype(str).str.strip())
 good=[]; rec=[]
 metal={3,4,11,12,13,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,87,88,89,90,91,92,93,94}
 for row in source.to_dict(orient='records'):
  uid=str(row.get('uid','')); rawcas=str(row.get('CAS No.','')).strip(); sm=row.get('Canonical_Smiles');remove=row.get('remove');outlier=row.get('outlier');
  reasons=[]
  if pd.notna(remove) and str(remove).strip():reasons.append('SOURCE_REMOVE_FLAG')
  if pd.notna(outlier) and str(outlier).strip():reasons.append('SOURCE_OUTLIER_FLAG')
  if not cas_valid(rawcas):reasons.append('CAS_INVALID_OR_NAME')
  lc=pd.to_numeric(row.get('Toxicity (mg/L)'),errors='coerce')
  if not np.isfinite(lc) or lc<=0:reasons.append('LC50_INVALID')
  m=Chem.MolFromSmiles(sm) if isinstance(sm,str) and sm.strip() else None
  if m is None: reasons.append('SMILES_INVALID')
  else:
   if len(Chem.GetMolFrags(m))!=1:reasons.append('MULTIFRAGMENT')
   if any(a.GetAtomicNum() in metal for a in m.GetAtoms()):reasons.append('METAL_HEAVY')
  can=Chem.MolToSmiles(m,isomericSmiles=True) if m else ''
  mw=float(Descriptors.MolWt(m)) if m else float('nan')
  if not np.isfinite(mw) or mw<=0:reasons.append('MW_INVALID')
  y=float(-np.log10(lc*0.001/mw)) if np.isfinite(lc) and lc>0 and np.isfinite(mw) and mw>0 else float('nan')
  baseline_ok=not any(x in reasons for x in ['SOURCE_REMOVE_FLAG','SOURCE_OUTLIER_FLAG','LC50_INVALID','SMILES_INVALID','MULTIFRAGMENT','METAL_HEAVY','MW_INVALID'])
  strict_ok=baseline_ok and cas_valid(rawcas)
  r={'uid':uid,'CAS':rawcas,'SMILES_original':row.get('SMILES'),'Canonical_Isomeric_SMILES':can,'LC50_mg_L':float(lc) if np.isfinite(lc) else np.nan,'MolWt_g_mol':mw,'pLC50_mol_L':y,'Original_pLD50_mgL_field':row.get('pLD50-final'),'Source_Remove_Flag':remove,'Source_Outlier_Flag':outlier,'CAS_Checkdigit_Pass':cas_valid(rawcas),'Baseline_Structure_Screen':bool(baseline_ok),'Strict_CAS_Screen':bool(strict_ok),'Reasons':';'.join(reasons),'In_Frozen_ECRA_CAS':rawcas in frozen_cas,'Original_Primary_Source_Verified':False,'Independent_Human_Review':False,'ECRA_QSAR_Released':False}
  rec.append(r)
 audit=pd.DataFrame(rec);audit.to_csv(OUT/'01_SOURCE_ROW_QC_1732.csv',index=False)
 base=audit[audit.Baseline_Structure_Screen].copy()
 strict=audit[audit.Strict_CAS_Screen].copy()
 # Deduplicate by canonical isomeric chemical identity, and collapse repeated labels using median.
 def collapse(x):
  vals=[]
  for sm,rows in x.groupby('Canonical_Isomeric_SMILES',sort=True):
   rr=rows.iloc[0]
   m=Chem.MolFromSmiles(sm)
   labels=rows.pLC50_mol_L.dropna().to_numpy()
   vals.append({'Chemical_ID':hashlib.sha1(sm.encode()).hexdigest()[:12],
     'Canonical_Isomeric_SMILES':sm,'CAS':';'.join(sorted(set(rows.CAS.astype(str)))),
     'Source_UIDs':';'.join(rows.uid.astype(str)),
     'n_source_rows':len(rows), 'y_range_log':float(labels.max()-labels.min()) if len(labels) else np.nan,
     'pLC50_mol_L':float(np.median(labels)),'LC50_mg_L_median':float(np.median(rows.LC50_mg_L)),
     'MolWt_g_mol':float(Descriptors.MolWt(m)),'Scaffold_Group':scaffold_of(m,sm),
     'Overlap_ECRA_CAS':bool(rows.In_Frozen_ECRA_CAS.any()),
     'Source_Study_ID_Available':False,'Approved_Training_Record':False})
  return pd.DataFrame(vals)
 bench=collapse(strict);bench.to_csv(OUT/'02_STRICT_CAS_CHEMICAL_BENCHMARK.csv',index=False)
 allchem=collapse(base);allchem.to_csv(OUT/'03_STRUCTURE_ONLY_SENSITIVITY_DATA.csv',index=False)
 assert len(bench)>=100, 'Insufficient strict CAS cohort'
 X,D,fps=split_features(bench);y=bench.pLC50_mol_L.to_numpy(dtype=float); groups=bench.Scaffold_Group.to_numpy(dtype=str)
 splitter=GroupShuffleSplit(n_splits=1,test_size=0.2,random_state=20261122)
 tr,te=next(splitter.split(X,y,groups=groups))
 assert set(groups[tr]).isdisjoint(set(groups[te]))
 assert set(bench.Canonical_Isomeric_SMILES.iloc[tr]).isdisjoint(set(bench.Canonical_Isomeric_SMILES.iloc[te]))
 splits=bench[['Chemical_ID','CAS','Canonical_Isomeric_SMILES','Scaffold_Group','pLC50_mol_L']].copy()
 splits['Partition']='TRAIN';splits.loc[te,'Partition']='SCAFFOLD_HELD_OUT';splits.to_csv(OUT/'04_LOCKED_SCAFFOLD_SPLIT.csv',index=False)
 Xtr=X[tr];Xte=X[te];ytr=y[tr];yte=y[te];gtr=groups[tr]
 cv=GroupKFold(n_splits=4)
 modeldefs=get_models();cvresults=[];cvpred=[]
 for name,model in modeldefs.items():
  for fold,(itr,iva) in enumerate(cv.split(Xtr,ytr,groups=gtr),start=1):
   m=model
   m.fit(Xtr[itr],ytr[itr]);pred=m.predict(Xtr[iva]);mm=metric(ytr[iva],pred)
   cvresults.append({'Model':name,'Fold':fold,'Test_n':len(iva),**mm})
   for ix,pr in zip(tr[iva],pred):cvpred.append({'Model':name,'Fold':fold,'Chemical_ID':bench.iloc[ix].Chemical_ID,'y_true':float(y[ix]),'y_pred':float(pr)})
 cvdf=pd.DataFrame(cvresults);cvdf.to_csv(OUT/'05_SCaffold_group_CV_metrics.csv',index=False)
 pd.DataFrame(cvpred).to_csv(OUT/'06_OUT_OF_FOLD_PREDICTIONS.csv',index=False)
 modelorder=[]
 for name,rows in pd.DataFrame(cvpred).groupby('Model'):
  modelorder.append((name,metric(rows.y_true.to_numpy(),rows.y_pred.to_numpy())))
 chosen=min(modelorder,key=lambda x:x[1]['RMSE'])[0]
 model=modeldefs[chosen];model.fit(Xtr,ytr)
 yp=model.predict(Xte);predtrain=model.predict(Xtr)
 joblib.dump({'estimator':model,'fingerprint_type':'Morgan radius 2 512 bit','descriptor_schema':['MolWt','LogP','TPSA','HBD','HBA','RotB','Rings','FractionCSP3','HeavyAtoms','AromaticRings','Heteroatoms','FormalCharge'],'source_flag':'EXPLORATORY_NOT_PRIMARY_EVIDENCE_APPROVED','selected_model':chosen,'model_training_dataset_sha256':sha(OUT/'02_STRICT_CAS_CHEMICAL_BENCHMARK.csv')},MOD/'EXPLORATORY_QSAR_NOT_APPROVED.joblib')
 # chemical nearest-neighbour applicability domain, train-defined threshold.
 trainfp=[fps[i] for i in tr];testfp=[fps[i] for i in te];nntrain=[]
 for i,f in enumerate(trainfp):
  vals=DataStructs.BulkTanimotoSimilarity(f,trainfp)
  vals[i]=-1
  nntrain.append(max(vals))
 threshold=float(np.quantile(nntrain,.05))
 nntest=[float(max(DataStructs.BulkTanimotoSimilarity(fp,trainfp))) for fp in testfp]
 tab=[]
 for k,ix in enumerate(te):
  rr=bench.iloc[ix]
  tab.append({'Chemical_ID':rr.Chemical_ID,'CAS':rr.CAS,'Canonical_Isomeric_SMILES':rr.Canonical_Isomeric_SMILES,'Scaffold_Group':rr.Scaffold_Group,'pLC50_Observed':float(yte[k]),'pLC50_Predicted':float(yp[k]),'Residual':float(yte[k]-yp[k]),'Max_Training_Tanimoto':nntest[k],'Inside_Training_Defined_AD':nntest[k]>=threshold})
 testdf=pd.DataFrame(tab);testdf.to_csv(OUT/'07_SCAFFOLD_HOLDOUT_TEST_PREDICTIONS.csv',index=False)
 # Group-bootstrap CIs for MAE/RMSE/R2 (held out scaffolds only).
 rng=np.random.default_rng(20261122);unig=np.unique(groups[te]);ci=[]
 if len(unig)>3:
  for _ in range(350):
   selected=rng.choice(unig,size=len(unig),replace=True)
   ids=np.concatenate([np.where(groups[te]==z)[0] for z in selected])
   if len(ids)>=3:ci.append(metric(yte[ids],yp[ids]))
 ci_summary={k:[float(np.nanpercentile([v[k] for v in ci if v[k] is not None],2.5)),float(np.nanpercentile([v[k] for v in ci if v[k] is not None],97.5))] for k in ('MAE','RMSE','R2')} if ci else {}
 # Train-only y scramble: fit model for fixed chosen estimator to shuffled labels and score on unchanged holdout; exploratory sanity check.
 perm=[]
 for i in range(20):
  shuffle=rng.permutation(ytr)
  # fast ridge for response permutation; use fixed not selected model
  mm=get_models()['Ridge'];mm.fit(Xtr,shuffle);pr=mm.predict(Xte)
  perm.append({'Permutation':i+1,**metric(yte,pr)})
 pd.DataFrame(perm).to_csv(OUT/'08_RESPONSE_PERMUTATION_SANITY_CHECK.csv',index=False)
 # Sensitivity: all parseable and non-flagged structures incl unreliable CAS; strict subset is primary.
 XX,_,_=split_features(allchem);yy=allchem.pLC50_mol_L.to_numpy(); gg=allchem.Scaffold_Group.to_numpy()
 aa,bb=next(GroupShuffleSplit(n_splits=1,test_size=.2,random_state=20261122).split(XX,yy,groups=gg))
 alt=get_models()[chosen];alt.fit(XX[aa],yy[aa]);yp_alt=alt.predict(XX[bb]);alt_perf=metric(yy[bb],yp_alt)
 sens={'Full_Structure_Cohort':len(allchem),'Strict_CAS_Cohort':len(bench),'All_Structures_Holdout':alt_perf}
 # Random chemical split sensitivity; no identical structures repeat, but scaffold overlap allowed.
 trr,ter=train_test_split(np.arange(len(bench)),test_size=.2,random_state=20261122)
 rr=get_models()[chosen];rr.fit(X[trr],y[trr]);p_rand=rr.predict(X[ter]);rand_perf=metric(y[ter],p_rand)
 # Distribution and diagnostics.
 metrics={'Dataset':'EXTERNAL SECONDARY BENCHMARK ONLY — NOT ECRA CONTROLLED QSAR',
 'Evidence_Caveat':'The source lacks per-row original experimental study identifiers and independent human authentication; no ECRA QSAR approval is asserted.',
 'External_Source':'https://github.com/zhaoqi106/ATFPGT-multi/blob/main/data_origin.csv',
 'External_Source_Git_Blob_SHA':'42a433d1722e0c0e7b904eeb3d618466458e8917',
 'Input_SHA256':sha(PINPUT),'Frozen_ECRA_SHA256':sha(FROZEN),'Row_Count_Source_RT':len(source),
 'Source_Excluded_or_Outlier':int((~audit['Baseline_Structure_Screen']).sum()),
 'Baseline_Chemical_Rows':len(base),'Invalid_CAS_After_Chemical_QC':int((base.Strict_CAS_Screen==False).sum()),
 'Strict_Row_Count':len(strict),'Strict_Chemical_Count':len(bench),'All_Valid_Structure_Chemical_Count':len(allchem),
 'Overlap_ECRA_CAS_Strict':int(bench.Overlap_ECRA_CAS.sum()),'Missing_Primary_Source_References':True,
 'Chosen_Model':chosen,'Feature_Count':X.shape[1],
 'Scaffold_Groups_Total':int(len(set(groups))),'Scaffold_Groups_Train':int(len(set(groups[tr]))),'Scaffold_Groups_Test':int(len(set(groups[te]))),
 'Scaffold_Groups_Overlap':int(len(set(groups[tr])&set(groups[te]))),
 'Train_n':len(tr),'Holdout_n':len(te),
 'Holdout':metric(yte,yp),'Train_In_Sample':metric(ytr,predtrain),
 'Train_CV_Models':{n:m for n,m in modelorder},
 'Holdout_Group_Bootstrap_95CI':ci_summary,
 'AD_Threshold_Training_NN_5thPercentile':threshold,
 'Holdout_AD_Coverage':float(np.mean(testdf.Inside_Training_Defined_AD)),
 'Holdout_AD_Inside_Metrics':metric(yte[testdf.Inside_Training_Defined_AD.to_numpy()],yp[testdf.Inside_Training_Defined_AD.to_numpy()]) if int(testdf.Inside_Training_Defined_AD.sum())>=3 else None,
 'Random_Chemical_Split_Sensitivity':rand_perf,
 'Structure_Only_Expanded_Sensitivity':sens,
 'Final_QSAR_Controlled_Release':'HOLD',
 'Unapproved_Predictive_Claims_Prohibited':True,
 'Software':{'rdkit':rdkit.__version__,'scikit_learn':sklearn.__version__},
 'Baseline_Median_Holdout':metric(yte,np.full_like(yte,np.median(ytr)))}
 (OUT/'09_QSAR_SUMMARY.json').write_text(json.dumps(metrics,indent=2),encoding='utf-8')
 # Charts
 plt.figure(figsize=(6.8,5.5))
 plt.scatter(yte,yp,alpha=.58,s=24)
 a=min(yte.min(),yp.min());b=max(yte.max(),yp.max());plt.plot([a,b],[a,b],linestyle='--',color='gray')
 plt.xlabel('Observed pLC50, molar (secondary dataset)');plt.ylabel('Predicted pLC50, molar');plt.title('Scaffold-held-out research benchmark (NOT validated ECRA QSAR)');plt.tight_layout();plt.savefig(FIG/'Fig76A_Scaffold_Holdout_Observed_vs_Predicted.png',dpi=250);plt.close()
 plt.figure(figsize=(7,5));plt.scatter(nntest,np.abs(yte-yp),alpha=.55,s=23)
 plt.axvline(threshold,linestyle='--',color='gray');plt.xlabel('Max. training-set Morgan Tanimoto');plt.ylabel('Absolute held-out error (pLC50)');plt.title('Research benchmark: similarity and prediction error');plt.tight_layout();plt.savefig(FIG/'Fig76B_Applicability_Domain.png',dpi=250);plt.close()
 # Verify frozen ECRA release remains untouched / no training original ECRA records
 assert frozen.Final_Eligibility.fillna('').astype(str).str.upper().eq('HOLD').all()
 assert not any(bench['Approved_Training_Record'])
 assert metrics['Scaffold_Groups_Overlap']==0
 assert not set(bench.Chemical_ID.iloc[tr]) & set(bench.Chemical_ID.iloc[te])
 print(json.dumps({k:metrics[k] for k in ['Row_Count_Source_RT','Strict_Chemical_Count','All_Valid_Structure_Chemical_Count','Train_n','Holdout_n','Chosen_Model','Holdout','Holdout_AD_Coverage','Random_Chemical_Split_Sensitivity','Final_QSAR_Controlled_Release']},indent=2))

if __name__=='__main__':main()
