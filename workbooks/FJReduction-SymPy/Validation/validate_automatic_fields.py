"""Exercise notebook field settings, including changes in the same kernel."""
from pathlib import Path
import os,sys,json,time
import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output
ROOT=Path(__file__).resolve().parents[1]
os.chdir(ROOT);sys.path.insert(0,str(ROOT))

def run_cell(shell,source):
    with capture_output() as captured:
        result=shell.run_cell(source,store_history=False)
    error=result.error_before_exec or result.error_in_exec
    if error: raise error
    return captured.stdout

records=[]
for name in ['MyExamples-Symbolic.ipynb','GettingStarted-Symbolic.ipynb']:
    nb=nbformat.read(ROOT/name,as_version=4)
    shell=InteractiveShell.instance();shell.reset(new_session=True)
    for cell in nb.cells:
        if cell.cell_type=='code': run_cell(shell,cell.source)
    for p,n in [(2,2),(3,2),(5,1)]:
        tick=time.perf_counter()
        run_cell(shell,f'p={p}\nn={n}\nmode="concrete"')
        for cell in nb.cells:
            if cell.cell_type=='code' and ('# Clear old results' in cell.source or 'if example_result is not None:' in cell.source or 'partial_comparison = None' in cell.source):
                run_cell(shell,cell.source)
        ns=shell.user_ns
        assert ns['example_result'] is not None,ns['calculation_error']
        assert ns['parameters'].q==p**n
        assert ns['sequence']['ExactSequenceVerified']
        records.append({'notebook':name,'p':p,'n':n,'q':p**n,'same_kernel_field_change':True,'seconds':round(time.perf_counter()-tick,3),'passed':True})
    run_cell(shell,'p=4\nn=1')
    computation=next(c.source for c in nb.cells if c.cell_type=='code' and '# Clear old results' in c.source)
    run_cell(shell,computation)
    assert shell.user_ns['example_result'] is None
    assert shell.user_ns['sequence'] is None
    assert shell.user_ns['calculation_error'] is not None
    run_cell(shell,'mode="symbolic"')
    run_cell(shell,computation)
    assert shell.user_ns['sequence']['ExactSequenceVerified']
    records.append({'notebook':name,'invalid_setting_clears_results':True,'formal_mode_passed':True})
for name in ['GettingStarted.ipynb','GettingStarted-Concrete.ipynb']:
    nb=nbformat.read(ROOT/name,as_version=4)
    shell=InteractiveShell.instance();shell.reset(new_session=True)
    run_cell(shell,nb.cells[1].source)
    run_cell(shell,nb.cells[3].source)
    for p,n in [(2,2),(3,2)]:
        code=nb.cells[3].source.replace('p = 3',f'p = {p}').replace('n = 1',f'n = {n}')
        run_cell(shell,code)
        for cell_index in [5,7,9,25]:
            run_cell(shell,nb.cells[cell_index].source)
        assert shell.user_ns['my_parameters'].q==p**n
        assert shell.user_ns['my_data'].parameters.q==p**n
        # Exact-sequence result returned by the same notebook's descriptive-report cell.
        run_cell(shell,nb.cells[29].source)
        assert shell.user_ns['walkthrough_result']['Sequence']['ExactSequenceVerified']
        records.append({'notebook':name,'p':p,'n':n,'q':p**n,'own_example_passed':True})
report={'version':'0.4.5','all_passed':True,'records':records}
(ROOT/'Validation/automatic-fields-0.4.5.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
