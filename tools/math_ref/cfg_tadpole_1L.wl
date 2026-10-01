(* 1-loop tadpole: tadpole[a] = ∫ d^d l / (l^2 - msq)^a.
   Closed form: I[1] = msq^(d/2 - 1) Γ(1 - d/2) = -Γ(-1+ε) * msq^(1-ε).
   We push epsorder high enough that truncation error at eps=0.01 is < 1e-15. *)
name        = "tadpole_1L";
family      = tad;
loops       = {l};
legs        = {};
conservation = {};
replacement = {};
propagators = {l^2 - msq};
numeric     = {msq -> 1};
targets     = {j[tad, 1], j[tad, 2]};
precision   = 50;
epsorder    = 12;
