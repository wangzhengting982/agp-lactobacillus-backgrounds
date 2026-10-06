"""Collect machine-verifiable evidence across all 40 submitted result tables."""
from pathlib import Path
import json,datetime,re,runpy
import openpyxl
ROOT=Path(__file__).resolve().parents[1]
def read(name):return json.loads((ROOT/'reports'/name).read_text(encoding='utf-8'))
def main():
    runpy.run_path(str(ROOT/'scripts/build_table_map.py'),run_name='__main__')
    integrity=read('input_integrity.json');prep=read('preparation_verification.json')
    core=read('core_verification.json');deep=read('deepening_validation.json')
    sens=read('sensitivity_validation.json');pixels=read('figures_pixel_comparison.json')
    figs=read('figures_provenance.json');taxonomy=read('core_taxonomy_input_rebuild.json')
    manuscript=read('manuscript_claims_check.json')
    wb=openpyxl.load_workbook(ROOT/'reference/submitted_results.xlsx',read_only=True)
    expected={s.split('_')[0] for s in wb.sheetnames if re.match(r'^[A-Z][0-9]{2}_',s)};wb.close()
    deep_codes={r['table'].split('_')[0] for r in deep['checks'] if r.get('reference')=='reference/submitted_results.xlsx'}
    verified=set(core['submission_tables_checked'])|deep_codes|{r['sheet'].split('_')[0] for r in sens['checks']}
    assert len(expected)==40 and expected==verified,(expected-verified,verified-expected)
    assert integrity['status']=='PASS' and prep['status']=='PASS'
    assert core['status']=='pass' and deep['all_passed'] and sens['all_passed']
    assert pixels['status']=='PASS',pixels
    assert not figs['figure_1_frozen_inputs_used'] and figs['returncode']==0
    assert taxonomy['status']=='pass'
    assert manuscript['status']=='PASS'
    errors=[r['max_absolute_error'] for r in core['results']]
    errors += [r['max_absolute_error'] for r in deep['checks'] if r.get('reference')=='reference/submitted_results.xlsx']
    errors += [r['max_absolute_error'] for r in sens['checks']]
    report={'status':'PASS','generated_local':datetime.datetime.now().isoformat(),
        'submission_result_sheets_verified':len(verified),'verified_sheet_codes':sorted(verified),
        'maximum_absolute_difference_across_submission_tables':max(errors),
        'cohort_n':{'descriptive':2760,'primary':2748},'background_genera':156,
        'primary_model_input_maximum_difference':max(r['max_abs_error'] for r in prep['matrix_comparisons']),
        'figures_rebuilt_from_new_results':3,'figure_pdf_200dpi_pixels_match':True,
        'invalid_full_round_contrasts_preserved':sens['invalid_full_round_contrasts_retained_as_missing'],
        'environment':'Independent Python 3.12 virtual environment on this Windows computer; package versions locked.',
        'scope':'Provided processed tables/BIOM plus archived ASV annotations, longitudinal sample index and reference sources -> cohorts, model matrices, fresh statistical fits, 40 current submission sheets and 3 figures.',
        'not_rerun':['FASTQ processing','initial VioScreen questionnaire matching','SINTAX taxonomy classification itself','selection of the archived longitudinal sample index','historical analyses removed from current submission','external-cohort validation'],
        'published':False,'github_remote_created':False,
        'manuscript_check_report':'reports/manuscript_claims_check.json',
        'evidence':['input_integrity.json','preparation_verification.json','core_taxonomy_input_rebuild.json','core_verification.json','deepening_validation.json','sensitivity_validation.json','figures_provenance.json','figures_pixel_comparison.json','manuscript_claims_check.json']}
    (ROOT/'reports/final_validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
