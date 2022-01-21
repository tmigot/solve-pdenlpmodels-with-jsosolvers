<!--This file was generated, do not modify it.-->
# Solve a PDE-constrained optimization problem with JSO-compliant solvers

In this tutorial you will learn how to use JSO-compliant solvers to solve a PDE-constrained optimization problem discretized with [PDENLPModels.jl](https://github.com/JuliaSmoothOptimizers/PDENLPModels.jl).

\toc

## Problem Statement

In this first part, we define a distributed Poisson control problem  with Dirichlet boundary conditions which is then automatically discretized.
We refer to [Gridap.jl](https://github.com/gridap/Gridap.jl) for more details on modeling PDEs and [PDENLPModels.jl](https://github.com/JuliaSmoothOptimizers/PDENLPModels.jl) for PDE-constrained optimization problems.

Let Ω = (-1,1)², we solve the following problem:

\begin{aligned}
\text{min}_{y \in H^1_0, u \in H^1} \quad &  \frac{1}{2} \int_\Omega |y(x) - yd(x)|^2dx + \frac{\alpha}{2} \int_\Omega |u|^2dx \\
\text{s.t.} & -\Delta y = h + u, \quad x \in \Omega, \\
& y = 0, \quad x \in \partial \Omega,
\end{aligned}

where yd(x) = -x₁² and α = 1e-2.
The force term is h(x₁, x₂) = - sin(ω x₁)sin(ω x₂) with  ω = π - 1/8.

To solve an optimization problem with Ipopt, the first thing to do is define your problem.
In this example, let's assume we want to solve the following problem:
\begin{aligned}
\text{min}_x \quad & (x_1 - 1)^2 + 4 (x_2 - x_1)^2 \\
\text{s.to} \quad & x_1^2 + x_2^2 \leq 1 \\
& x_1 \leq 0.5 \\
& 0.25 \leq x_2 \leq 0.75
\end{aligned}

```julia:ex1
using Gridap, PDENLPModels
```

Definition of the domain and discretization

```julia:ex2
n = 3
domain = (-1, 1, -1, 1)
partition = (n, n)
model = CartesianDiscreteModel(domain, partition)
```

Definition of the FE-spaces

```julia:ex3
reffe = ReferenceFE(lagrangian, Float64, 2)
Xpde = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
y0(x) = 0.0
Ypde = TrialFESpace(Xpde, y0)

reffe_con = ReferenceFE(lagrangian, Float64, 1)
Xcon = TestFESpace(model, reffe_con; conformity = :H1)
Ycon = TrialFESpace(Xcon)
Y = MultiFieldFESpace([Ypde, Ycon])
```

Integration machinery

```julia:ex4
trian = Triangulation(model)
degree = 1
dΩ = Measure(trian, degree)
```

Objective function

```julia:ex5
yd(x) = -x[1]^2
α = 1e-2
function f(y, u)
  ∫(0.5 * (yd - y) * (yd - y) + 0.5 * α * u * u) * dΩ
end
```

Definition of the constraint operator

```julia:ex6
ω = π - 1 / 8
h(x) = -sin(ω * x[1]) * sin(ω * x[2])
function res(y, u, v)
  ∫(∇(v) ⊙ ∇(y) - v * u - v * h) * dΩ
end
op = FEOperator(res, Y, Xpde)
```

Definition of the initial guess

```julia:ex7
npde = Gridap.FESpaces.num_free_dofs(Ypde)
ncon = Gridap.FESpaces.num_free_dofs(Ycon)
x0 = zeros(npde + ncon);
```

Overall, we built a GridapPDENLPModel, which implements the [NLPModel](https://juliasmoothoptimizers.github.io/NLPModels.jl/stable/) API.

```julia:ex8
nlp = GridapPDENLPModel(x0, f, trian, Ypde, Ycon, Xpde, Xcon, op, name = "Control elastic membrane")

using NLPModels

(get_nvar(nlp), get_ncon(nlp))
```

# ## Find a Feasible Point

Before solving the previously defined model, we will first improve our initial guess.
The first step is to create a nonlinear least-squares whose residual is the equality-constraint of the optimization problem.
We use `FeasibilityResidual` from [NLPModelsModifiers.jl](https://github.com/JuliaSmoothOptimizers/NLPModelsModifiers.jl) to convert the NLPModel as an NLSModel.
Then, using `trunk`, a matrix-free solver for least-squares problems implemented in [JSOSolvers.jl](https://github.com/JuliaSmoothOptimizers/JSOSolvers.jl), we find an
improved guess which is close to being feasible for our large-scale problem.
By default, JSO-compliant solvers use `nlp.meta.x0` as an initial guess.

```julia:ex9
using JSOSolvers, NLPModelsModifiers

nls = FeasibilityResidual(nlp)
stats_trunk = trunk(nls)
```

We check the solution from the stats returned by `trunk`:

```julia:ex10
norm(cons(nlp, stats_trunk.solution))
```

We will use the solution found to initialize our solvers.

## Solve the Problem

Finally, we are ready to solve the PDE-constrained optimization problem with a targeted tolerance of `10⁻⁵`.
In the following, we will use both Ipopt and DCI on our problem.

```julia:ex11
using NLPModelsIpopt
```

Set `print_level = 0` to avoid printing detailed iteration information.

```julia:ex12
stats_ipopt = ipopt(nlp, x0 = stats_trunk.solution, tol = 1e-5, print_level = 0)
```

The problem was successfully solved, and we can extract the function evaluations from the stats.

```julia:ex13
stats_ipopt.counters
```

Reinitialize the counters before re-solving.

```julia:ex14
reset!(nlp);
```

Most JSO-compliant solvers are using logger for printing iteration information.
`NullLogger` avoids printing iteration information.

```julia:ex15
using DCISolver, Logging

stats_dci = with_logger(NullLogger()) do
  dci(nlp, stats_trunk.solution, atol = 1e-5, rtol = 0.0)
end
```

The problem was successfully solved, and we can extract the function evaluations from the stats.

```julia:ex16
stats_dci.counters
```

We now compare the two solvers with respect to the time spent,

```julia:ex17
stats_ipopt.elapsed_time, stats_dci.elapsed_time
```

and also check objective value, feasibility, and dual feasibility of `ipopt` and `dci`.

```julia:ex18
(stats_ipopt.objective, stats_ipopt.primal_feas, stats_ipopt.dual_feas),
(stats_dci.objective, stats_dci.primal_feas, stats_dci.dual_feas)
```

Overall `DCISolver` is doing great for solving large-scale optimization problems!
You can try increase the problem size by changing the discretization parameter `n`.

