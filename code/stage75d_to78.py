#!/usr/bin/env python3
"""ECRA Stage75D-78 conservative research release and feasibility runner.

Never signs evidence or releases model records without actual independent, recorded
reviewer and adjudicator decisions.  The unmodified frozen input remains separate.
"""
import csv, json, hashlib, collections, os, math, datetime
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
I=ROOT/'input'; O=ROOT/'output'; O.mkdir(exist_ok=True)
VALID={'APPROVE','REJECT','HOLD'}
FROZEN=I/'frozen_169.csv'
EXCLUDED_SPECIES={'37424','37960'}
FORMULATIONS={'100091','116077'}
TRANSFORMATION={'37959','37423'}
TAUTOMER={'125793'}
SOURCE_PDFS={'13142':'13142_GALASSI_1988.pdf','569':'569_DEGRAEVE_1980.pdf','12605':'12605_WAN_1987.pdf'}
PRIORITY={'13142':1,'569':2,'12605':3}


def read(path):
 with open(path,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(path,rows,fields=None):
 if not fields:
  fields=list(rows[0]) if rows else []
 with open(path,'w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for chunk in iter(lambda:f.read(65536),b''):h.update(chunk)
 return h.hexdigest()
def dateok(s):
 try:
  d=datetime.datetime.fromisoformat(s.replace('Z','+00:00'))
  return d.tzinfo is not None
 except Exception:return False

def template_rows():
 frozen={r['Result_ID']:r for r in read(FROZEN)}
 chem={r['CAS']:r for r in read(I/'stage59_structure_107.csv')}
 st44={r['Result_ID']:r for r in read(I/'stage44_descriptors.csv')}
 prim=read(I/'primary_22.csv');out=[]
 for p in sorted(prim,key=lambda z:(PRIORITY.get(z['Reference_Number'],9),z['Result_ID'])):
  q=frozen[p['Result_ID']]
  reasons=[]
  if p['Target_Species_Compatible']!='YES':reasons.append('WRONG_SPECIES_HARD_EXCLUSION')
  if p['Result_ID'] in FORMULATIONS:reasons.append('COMMERCIAL_FORMULATION_NOT_SINGLE_MOLECULE')
  if p['Result_ID'] in TRANSFORMATION:reasons.append('ACUTE_TEST_TRANSFORMATION_OR_STABILITY_REVIEW')
  if p['Result_ID'] in TAUTOMER:reasons.append('TAUTOMER_NORMALIZATION_REVIEW')
  if p['CAS'] in {'71-43-2','108-88-3','95-47-6'}:reasons.append('OTHER_REFERENCE_CAS_OVERLAP_UNVERIFIED')
  sr=st44.get(p['Result_ID'],{})
  r={'Result_ID':p['Result_ID'],'Test_ID':p['Test_ID'],'Reference_Number':p['Reference_Number'],
    'CAS':p['CAS'],'Original_Test_Material':p['Original_Test_Material'],
    'Primary_Table_LC50_mg_L':p['Original_LC50_mg_L'],'Original_Page':p['Study_Page_Printed'],
    'Study_DOI':p['Study_DOI'],'Original_Source_Species':p['Original_Source_Species'],
    'Original_Table':p['Study_Table'], 'Source_Table_Numerical_Match':p['Numerical_Table_Match'],
    'Target_Trout_Compatible':p['Target_Species_Compatible'],
    'Stage59_Screen':chem.get(p['CAS'],{}).get('Structure_Screen','UNKNOWN'),
    'Stage44_SMILES_PROVISIONAL':sr.get('Canonical_SMILES',''),
    'Stage44_MW_PROVISIONAL':sr.get('MolWt',''),
    'Specific_Risk_Flags':'|'.join(reasons) or 'NO_SPECIFIC_RISK_DETECTED_NOT_APPROVED',
    'Review_Priority':'CRITICAL' if reasons else 'STANDARD',
    'Current_Release':'HOLD',
    'Verified_Source_PDF_SHA256':'',
    'Evidence_Page_Reference':'',
    'Reviewed_Test_Form_Description':'',
    'Reviewed_Analyte_Concentration_Basis':'',
    'Reviewed_Canonical_Isomeric_SMILES':'',
    'Structure_Change_Justification':'', 'Other_Reference_Lineage_Resolution':'',
    'Reviewed_Molecular_Weight_g_mol':'',
    'Source_Lineage_Independence_Explanation':'',
    'Transformation_Resolution':'',
    'Reviewer_A_Name':'','Reviewer_A_UTC':'','Reviewer_A_Decision':'',
    'Reviewer_A_Evidence_Justification':'','Reviewer_B_Name':'','Reviewer_B_UTC':'',
    'Reviewer_B_Decision':'','Reviewer_B_Evidence_Justification':'',
    'Identity_Gate':'','Analyte_Basis_Gate':'','Study_Independence_Gate':'',
    'Final_Adjudicator_Name':'','Final_Adjudication_UTC':'',
    'Final_Adjudication_Decision':'','Final_Adjudicator_Justification':''}
  out.append(r)
 return out

def evaluate(row):
 rid=row['Result_ID'];errors=[]
 if rid in EXCLUDED_SPECIES: return 'EXCLUDE_WRONG_SPECIES',['SOURCE_TABLE_SPECIES_MISMATCH']
 if rid in FORMULATIONS: return 'EXCLUDE_MULTICOMPONENT_FORMULATION',['NOT_SINGLE_MOLECULE_EXPOSURE']
 if row['Source_Table_Numerical_Match']!='YES' or row['Target_Trout_Compatible']!='YES':errors.append('PRIMARY_SOURCE_TABLE_AND_SPECIES_NOT_AUTHENTICATED')
 if row['Stage59_Screen']!='COMPUTATIONALLY_ELIGIBLE_PROVISIONAL':errors.append('NO_PARSEABLE_SINGLE_STRUCTURE')
 if not row['Verified_Source_PDF_SHA256'].strip(): errors.append('PDF_HASH_NOT_ATTESTED')
 if not row['Evidence_Page_Reference'].strip():errors.append('TABLE_PAGE_NOT_ATTESTED')
 for k in ('Reviewed_Test_Form_Description','Reviewed_Analyte_Concentration_Basis','Reviewed_Canonical_Isomeric_SMILES','Reviewed_Molecular_Weight_g_mol','Source_Lineage_Independence_Explanation'):
  if not row[k].strip(): errors.append(k+'_MISSING')
 try:
  mw=float(row['Reviewed_Molecular_Weight_g_mol'])
  if not (20<mw<2000 and math.isfinite(mw)):errors.append('IMPLAUSIBLE_OR_INVALID_MW')
 except Exception: errors.append('UNVERIFIED_MW')
 # A parseable SMILES and a form-specific molecular-weight match are required;
 # provisional Stage44 structure alone is never an approval.
 try:
  from rdkit import Chem
  from rdkit.Chem import Descriptors
  smi=row['Reviewed_Canonical_Isomeric_SMILES'].strip()
  mol=Chem.MolFromSmiles(smi)
  if mol is None:errors.append('INVALID_REVIEWED_SMILES')
  else:
   if not math.isclose(Descriptors.MolWt(mol),float(row['Reviewed_Molecular_Weight_g_mol']),rel_tol=0.0,abs_tol=0.2):errors.append('VERIFIED_SMILES_MW_DISAGREEMENT')
   old=Chem.MolFromSmiles(row['Stage44_SMILES_PROVISIONAL'].strip())
   if old is not None and Chem.MolToSmiles(old,isomericSmiles=False)!=Chem.MolToSmiles(mol,isomericSmiles=False):
    if len(row.get('Structure_Change_Justification','').strip())<60:errors.append('STRUCTURE_CHANGE_REQUIRES_SOURCE_EXPLANATION')
 except ImportError:errors.append('RDKIT_REQUIRED_FOR_MODEL_RELEASE')
 except Exception:errors.append('CHEMICAL_IDENTITY_VALIDATION_FAILED')
 # The 169-observation register is the authority for identifying CAS re-use.
 cas_references={x['Reference_Number'] for x in read(FROZEN) if x['CAS']==row['CAS']}
 if len(cas_references)>1 and len(row.get('Other_Reference_Lineage_Resolution','').strip())<60:
  errors.append('OTHER_REFERENCE_SAME_CAS_LINEAGE_UNRESOLVED')
 for k in ('Identity_Gate','Analyte_Basis_Gate','Study_Independence_Gate'):
  if row[k]!='APPROVE':errors.append(k+'_NOT_APPROVED')
 if rid in TRANSFORMATION:
  if len(row['Transformation_Resolution'].strip())<60:errors.append('TRANSFORMATION_NOT_EVIDENCE_RESOLVED')
 if rid in TAUTOMER and len(row['Reviewed_Test_Form_Description'].strip())<30:errors.append('TAUTOMER_NOT_EXPLAINED')
 for i in 'AB':
  if row[f'Reviewer_{i}_Decision']!='APPROVE':errors.append('REVIEWER_'+i+'_NOT_APPROVED')
  if not row[f'Reviewer_{i}_Name'].strip() or not dateok(row[f'Reviewer_{i}_UTC']):errors.append('REVIEWER_'+i+'_IDENTITY_OR_UTC_MISSING')
  if len(row[f'Reviewer_{i}_Evidence_Justification'].strip())<60:errors.append('REVIEWER_'+i+'_RATIONALE_INSUFFICIENT')
 for n in ('Reviewer_A_Name','Reviewer_B_Name','Final_Adjudicator_Name'):
  if not row[n].strip():errors.append(n+'_MISSING')
 people=[row[k].strip().casefold() for k in ('Reviewer_A_Name','Reviewer_B_Name','Final_Adjudicator_Name')]
 if '' not in people and len(set(people))!=3:errors.append('REVIEWERS_AND_ADJUDICATOR_NOT_DISTINCT')
 if row['Final_Adjudication_Decision']!='APPROVE' or not dateok(row['Final_Adjudication_UTC']):errors.append('FINAL_ADJUDICATOR_NOT_SIGNED')
 if len(row['Final_Adjudicator_Justification'].strip())<60:errors.append('FINAL_RATIONALE_INSUFFICIENT')
 if errors:return 'HOLD_PENDING_INDEPENDENT_ADJUDICATION',sorted(set(errors))
 return 'CONDITIONAL_HUMAN_ATTESTED_RELEASE_FOR_FURTHER_QA',[]

def run():
 frozen=read(FROZEN); primary=read(I/'primary_22.csv')
 assert len(frozen)==169 and len(primary)==22
 assert len({r['Result_ID'] for r in frozen})==169
 form_path=O/'01_DUAL_REVIEW_AND_RELEASE_TEMPLATE_22.csv'
 # Preserve real user reviewer responses if present: never overwrite review decisions.
 if form_path.exists():
  form=read(form_path)
  for row in form:
   row.setdefault('Structure_Change_Justification','')
   row.setdefault('Other_Reference_Lineage_Resolution','')
  write(form_path,form)
  if set(r['Result_ID'] for r in form)!=set(r['Result_ID'] for r in primary):raise ValueError('Reviewer file no longer matches 22 primary records; manual reconcile required')
 else:
  form=template_rows();write(form_path,form)
 reviewed=[]
 for row in form:
  status,problems=evaluate(row)
  reviewed.append({**row,'Machine_Gate_Decision':status,'Blocking_Reasons':'|'.join(problems)})
 write(O/'02_22_RECORD_AUTOMATED_RELEASE_GATE.csv',reviewed)
 by_id={r['Result_ID']:r for r in reviewed};ready=[];gate=[]
 for f in frozen:
  rid=f['Result_ID'];e=by_id.get(rid)
  status=e['Machine_Gate_Decision'] if e else 'HOLD_NO_PRIMARY_TABLE'
  out={k:f[k] for k in ('Result_ID','Test_ID','Reference_Number','CAS','Species','Endpoint','Duration','Candidate_LC50_mg_L','Provisional_pLC50','Final_Eligibility','Model_Use')}
  out['Stage75D_Release_Gate']=status
  out['Stage75D_Approved']='YES' if status.startswith('CONDITIONAL_HUMAN') else 'NO'
  out['Stage75D_Blockers']=e['Blocking_Reasons'] if e else 'PRIMARY_TABLE_NOT_AUTHENTICATED'
  gate.append(out)
  if out['Stage75D_Approved']=='YES':ready.append(out)
 write(O/'03_ALL169_CONTROLLED_READONLY_GATE.csv',gate)
 fields=['Result_ID','Test_ID','Reference_Number','CAS','Species','Endpoint','Duration','Candidate_LC50_mg_L','Approved_Molecular_Weight_g_mol','Approved_SMILES','Newly_Calculated_pLC50','Approval_Provenance']
 approved=[]
 for r in reviewed:
  if r['Machine_Gate_Decision'].startswith('CONDITIONAL_HUMAN'):
   mg=float(r['Primary_Table_LC50_mg_L']);mw=float(r['Reviewed_Molecular_Weight_g_mol']);pl=-math.log10(mg*1e-3/mw)
   approved.append({'Result_ID':r['Result_ID'],'Test_ID':r['Test_ID'],'Reference_Number':r['Reference_Number'],'CAS':r['CAS'],'Species':'Oncorhynchus mykiss','Endpoint':'96h LC50','Duration':'96 h','Candidate_LC50_mg_L':mg,'Approved_Molecular_Weight_g_mol':mw,'Approved_SMILES':r['Reviewed_Canonical_Isomeric_SMILES'],'Newly_Calculated_pLC50':f'{pl:.9f}','Approval_Provenance':r['Final_Adjudicator_Name']})
 write(O/'04_CONDITIONALLY_RELEASED_DATASET.csv',approved,fields)
 valid_trout=[p for p in primary if p['Target_Species_Compatible']=='YES'];nonform=[p for p in valid_trout if p['Result_ID'] not in FORMULATIONS]
 stable=[p for p in nonform if p['Result_ID'] not in TRANSFORMATION]
 stats={'frozen_rows':len(frozen),'unique_CAS':len(set(q['CAS'] for q in frozen)),'references':len(set(q['Reference_Number'] for q in frozen)),
 'original_numerical_matches':len(primary),'trout_compatible_matches':len(valid_trout),'eligible_without_commercial_formulations_MAX_NOT_APPROVED':len(nonform),
 'without_formulations_or_two_reactive_compounds_MAX_NOT_APPROVED':len(stable),
 'source_groups_in_trout_compatible':len(set(p['Reference_Number'] for p in valid_trout)),
 'authoritative_reviewers_signed':sum(1 for r in reviewed if r['Machine_Gate_Decision'].startswith('CONDITIONAL_HUMAN')),
 'conditionally_released_rows':len(approved),'approved_for_QSAR':0,
 'validation_feasibility':'HOLD_NO_MODEL: no fully quality-certified source-independent training and external test',
 'project_policy_minimum_suggested_unique_chemicals':40,'project_policy_minimum_suggested_independent_sources':5,
 'policy_note':'Conservative prospective project-design thresholds only; not universally required by QSAR standards.',
 'data_read_only':True,'original_frozen_sha256':sha(FROZEN)}
 (O/'05_STAGE75D_78_FINAL_STATUS.json').write_text(json.dumps(stats,indent=2)+'\n',encoding='utf-8')
 actions=[]
 for r in reviewed:
  if r['Machine_Gate_Decision'].startswith('HOLD'):
   actions.append({'Reference':r['Reference_Number'],'Result_ID':r['Result_ID'],'CAS':r['CAS'],'Priority':r['Review_Priority'],'Decision':r['Machine_Gate_Decision'],'Action':'Source-form/analyte/lineage verification; two independent reviewer attestations; third adjudication','Reason':r['Blocking_Reasons']})
 write(O/'06_HUMAN_ACTION_REGISTER.csv',actions)
 hashes=[]
 for p in sorted(I.iterdir()):
  if p.is_file():hashes.append({'input':p.name,'sha256':sha(p),'bytes':p.stat().st_size})
 write(O/'07_INPUT_SHA256_MANIFEST.csv',hashes)
 print(json.dumps(stats,indent=2));return stats
if __name__=='__main__':run()
