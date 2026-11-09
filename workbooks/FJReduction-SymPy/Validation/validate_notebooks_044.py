"""Execute every code cell in a fresh IPython/Sage process, capture rich output,
and test changed inputs. Sage mode uses the Sage preparser. This verifies code
and emitted MIME output, not a Jupyter browser UI.
"""
from pathlib import Path
import argparse,contextlib,hashlib,json,os,subprocess,sys,time
import nbformat

KIND='sympy'
ROOT=Path(__file__).resolve().parents[1]

def execute(name,scenario=None):
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    os.chdir(ROOT);sys.path.insert(0,str(ROOT))
    if KIND=='sage':
        from sage.repl.interpreter import get_test_shell
        shell=get_test_shell()
    else:
        InteractiveShell.clear_instance();shell=InteractiveShell.instance()
    notebook=nbformat.read(ROOT/name,as_version=4)
    records=[];tick=time.perf_counter()
    for index,cell in enumerate(notebook.cells):
        if cell.cell_type!='code':continue
        code=cell.source
        if scenario and '# EDIT ONLY THESE TWO MATRICES' in code:
            code={
             'scalar':'Phi=[[theta+tau**2]]\nPsi=[[theta+tau]]',
             'nonscalar':'Phi=[[theta+tau**2,1],[0,theta+tau**2]]\nPsi=[[theta,1+tau],[0,theta]]',
             'new-coefficients':('Phi=[["theta+u*tau**3"]]\nPsi=[["theta+v*tau"]]' if KIND=='sympy' else 'Phi=[[theta+var("u")*tau**3]]\nPsi=[[theta+var("v")*tau]]'),
             'invalid':'Phi=[[theta+tau]]\nPsi=[[theta+tau]]'}[scenario]
        start=time.perf_counter()
        with capture_output() as captured:
            result=InteractiveShell.run_cell(shell,code,store_history=True)
        error=result.error_before_exec or result.error_in_exec
        if error:raise RuntimeError(f'{name}, cell {index+1}: {error}\n{captured.stdout}\n{captured.stderr}') from error
        cell.execution_count=len(records)+1;cell.outputs=[]
        if captured.stdout:cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
        if captured.stderr:cell.outputs.append(nbformat.v4.new_output('stream',name='stderr',text=captured.stderr))
        for output in captured.outputs:cell.outputs.append(nbformat.v4.new_output('display_data',data=output.data,metadata=output.metadata))
        records.append({'cell':index+1,'seconds':round(time.perf_counter()-start,4)})
    ns=shell.user_ns
    summary=ns['j2Counterexample']['Summary']
    assert tuple(summary[k] for k in ['FullDimension','TrueExt0Dimension','H','IndividualZeroCount','InvalidCountPrediction'])==(8,6,2,1,7)
    assert summary['PartialBasisStatus']=='Unavailable'
    assert ns['j2Counterexample']['Result']['Sequence']['ExactSequenceVerified']
    view=ns['j2Counterexample']['Output'];html=view._repr_html_()
    assert 'τ' in html and 'Why partial basis is unavailable' in html
    assert 'rank D' in view.text and 'Matrix size: 8 x 8' in view.text and 'Matrix size: 6 x 6' in view.text
    if name.startswith('MyExamples'):
        value=ns['example_result']
        if scenario=='invalid':assert value is None and ns['sequence'] is None and ns['calculation_error'] is not None
        else:
            expected={'scalar':(2,1,1),'nonscalar':(8,6,2),'new-coefficients':(3,2,1)}.get(scenario,(14,12,2))
            assert tuple(value['Sequence'][k] for k in ['FullDimension','ExtTau0Dimension','H'])==expected
            assert value['Sequence']['ExactSequenceVerified']
    if 'GettingStarted-Concrete' in name or name=='GettingStarted.ipynb' or (KIND=='sage' and 'GettingStarted-Symbolic' in name):
        path=ROOT/'Validation'/('native-validation.json' if KIND=='sage' else 'test-results.json')
        tests=json.loads(path.read_text());assert tests.get('successful',tests.get('all_passed'))
        assert tests.get('tests_run',tests.get('tests'))==(142 if KIND=='sage' else 120)
    # Deliver default notebooks with real executed output. Scenarios never replace input cells.
    nbformat.validate(notebook)
    if not scenario:nbformat.write(notebook,ROOT/name)
    result={'notebook':name,'scenario':scenario,'code_cells_executed':len(records),'all_passed':True,
      'seconds':round(time.perf_counter()-tick,3),'cell_timings':records,'rich_report_verified':True,
      'source_sha256':hashlib.sha256((ROOT/name).read_bytes()).hexdigest()}
    target=ROOT/'Validation'/(name.removesuffix('.ipynb')+'-'+(scenario or 'default')+'-execution.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(name,scenario or 'default',len(records),'cells passed',flush=True)
    return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--one');parser.add_argument('--scenario');args=parser.parse_args()
    if args.one:return execute(args.one,args.scenario)
    names=['GettingStarted-Concrete.ipynb','GettingStarted-Symbolic.ipynb','MyExamples-Symbolic.ipynb','PartialReduction-Example.ipynb','Counterexample-J2.ipynb']
    if KIND=='sympy':names.append('GettingStarted.ipynb')
    requests=[(n,None) for n in names]+[('MyExamples-Symbolic.ipynb',s) for s in ['scalar','nonscalar','new-coefficients','invalid']]
    records=[]
    for name,scenario in requests:
        cmd=[sys.executable,__file__,'--one',name]+(['--scenario',scenario] if scenario else [])
        print('Checking',name,scenario or 'default',flush=True)
        subprocess.run(cmd,cwd=ROOT,check=True)
        target=ROOT/'Validation'/(name.removesuffix('.ipynb')+'-'+(scenario or 'default')+'-execution.json')
        records.append(json.loads(target.read_text()))
    implementation='SymPy'
    runtime={}
    if KIND=='sage':
        from sage.version import version
        from importlib import metadata
        try:runtime['passagemath_package_version']=metadata.version('passagemath-standard');implementation='passagemath (compiled SageMath-compatible runtime)'
        except metadata.PackageNotFoundError:implementation='SageMath'
        runtime['sage_api_version']=str(version)
    result={'method':'Fresh processes with actual IPython code-cell execution and captured rich display',
      'implementation':implementation,**runtime,'jupyter_browser_ui_checked':False,
      'notebooks':records,'code_cells_executed':sum(r['code_cells_executed'] for r in records),'all_passed':True}
    (ROOT/'Validation/notebook-execution-0.4.5.json').write_text(json.dumps(result,indent=2)+'\n')
    print('All notebooks and input-change scenarios passed:',result['code_cells_executed'],'cells',flush=True)
    return result

if __name__=='__main__':main()
