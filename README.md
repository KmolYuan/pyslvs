[![PyPI](https://img.shields.io/pypi/v/pyslvs.svg)](https://pypi.org/project/pyslvs/)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/KmolYuan/pyslvs.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/KmolYuan/pyslvs/context:python)

# Pyslvs Libraries

A no-GUI module of mechanism synthesis system and a 2D geometric constraint solver.

## Installation

Install from PyPI:

```bash
pip install pyslvs
```

Or install and test from source:

```bash
pip install -e .
python test
```

## Documentation

Run the solver through an example:

```python
from pyslvs import example_list, parse_vpoints, t_config, expr_solving

# Get example with name
expr, inputs = example_list("Jansen's linkage (Single)")
# Parse the mechanism expression into a list of joint data
vpoints = parse_vpoints(expr)
# Config joint data and control data for the solver
exprs = t_config(vpoints, inputs)
# Solve the position
result = expr_solving(exprs, vpoints, {pair: 0. for pair in inputs})
# Get the result from joint 7
x, y = result[7]
print(x, y)  # -43.170055 -91.753226
```

The documentation of Pyslvs library is on [Readthedocs](https://pyslvs-ui.readthedocs.io/en/latest/pyslvs-lib/).

If you have any questions, please post on GitHub issue or contact <pyslvs@gmail.com>.
