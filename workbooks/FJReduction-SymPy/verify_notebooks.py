"""Run the walkthroughs in a fresh in-process IPython shell per notebook.

This records real code-cell execution and rich display, without requiring
Jupyter kernel sockets. It is not a Jupyter browser UI validation.
"""
from pathlib import Path
import hashlib
import json
import time
import nbformat
from IPython.core.interactiveshell import InteractiveShell
from IPython.utils.capture import capture_output

ROOT = Path(__file__).resolve().parent


def execute(name, overrides=None, write=True):
    notebook = nbformat.read(ROOT/name, as_version=4)
    InteractiveShell.clear_instance()
    shell = InteractiveShell.instance()
    start = time.perf_counter()
    count = 0
    for index, cell in enumerate(notebook.cells):
        if cell.cell_type != 'code':
            continue
        source = (overrides or {}).get(index, cell.source)
        with capture_output() as captured:
            result = shell.run_cell(source, store_history=True)
        error = result.error_before_exec or result.error_in_exec
        if error:
            raise RuntimeError(f'{name}: cell {index+1}: {error}\n{captured.stdout}\n{captured.stderr}') from error
        cell.execution_count = count+1
        cell.outputs = []
        if captured.stdout:
            cell.outputs.append(nbformat.v4.new_output('stream', name='stdout', text=captured.stdout))
        if captured.stderr:
            cell.outputs.append(nbformat.v4.new_output('stream', name='stderr', text=captured.stderr))
        for output in captured.outputs:
            cell.outputs.append(nbformat.v4.new_output('display_data', data=output.data,
                                                     metadata=output.metadata))
        count += 1
    nbformat.validate(notebook)
    if write:
        nbformat.write(notebook, ROOT/name)
    return {'notebook': name, 'code_cells': count, 'all_passed': True,
            'elapsed_seconds': round(time.perf_counter()-start, 3),
            'sha256': hashlib.sha256((ROOT/name).read_bytes()).hexdigest()}


def main():
    results = [execute('GettingStarted-Concrete.ipynb'), execute('MyExamples-Symbolic.ipynb'), execute('PartialReduction-Example.ipynb')]
    # A fresh run changes ONLY Phi/Psi to new coefficient names.
    n = nbformat.read(ROOT/'MyExamples-Symbolic.ipynb', as_version=4)
    input_index = next(i for i, cell in enumerate(n.cells)
                       if cell.cell_type == 'code' and '# EDIT ONLY THESE TWO MATRICES' in cell.source)
    results.append(execute('MyExamples-Symbolic.ipynb', {input_index:
                          'Phi = [["theta + u*tau**3"]]\nPsi = [["theta + v*tau"]]'}, write=False))
    results[-1]['scenario'] = 'Only Phi/Psi changed; undeclared symbols u,v discovered automatically'
    results.append(execute('MyExamples-Symbolic.ipynb', {input_index:
                          'Phi = [["theta + tau**2", 1], [0, "theta + tau**2"]]\nPsi = [["theta + tau", 1], [0, "theta + tau"]]'}, write=False))
    results[-1]['scenario'] = 'Only Phi/Psi changed; nonscalar quotient with matching J2 nilpotents'
    # Compatibility notebook is an exact copy of the verified concrete one.
    (ROOT/'GettingStarted.ipynb').write_bytes((ROOT/'GettingStarted-Concrete.ipynb').read_bytes())
    report = {'method': 'Fresh in-process IPython shells with captured rich display',
              'jupyter_browser_ui_checked': False, 'results': results,
              'total_executed_code_cells': sum(x['code_cells'] for x in results),
              'all_passed': all(x['all_passed'] for x in results)}
    (ROOT/'Validation/notebook-execution.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
