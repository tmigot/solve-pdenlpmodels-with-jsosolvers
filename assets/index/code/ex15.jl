# This file was generated, do not modify it. # hide
using DCISolver, Logging

stats_dci = with_logger(NullLogger()) do
  dci(nlp, stats_trunk.solution, atol = 1e-5, rtol = 0.0)
end