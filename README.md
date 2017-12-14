# sylie
Explicit SO(3) SE(3) and SIM(3) Jacobians implemented with sympy

A useful Jacobian for error reduction on lie manifolds is

```
d      
--  log(exp(a)^-1 * exp(b)) - error
db     
```

This repository is meant as an accompaniment to a blog post which
derives explicit formulae for SO(3) SE(3) and SIM(3) Jacobians of the above form.

