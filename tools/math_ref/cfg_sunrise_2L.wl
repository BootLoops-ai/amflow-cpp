(* 2-loop sunrise (3 equal-mass propagators).
   Tests: SingleMass real recursion (multi-loop vacuum). *)
name        = "sunrise_2L";
family      = sunrise;
loops       = {l1, l2};
legs        = {p};
conservation = {};
replacement = {p^2 -> s};
(* Pad with 2 auxiliary (irreducible) scalar products so the basis is complete.
   The physical sunrise master corresponds to indices {1,1,1,0,0}. *)
propagators = {l1^2 - msq, l2^2 - msq, (l1 + l2 - p)^2 - msq,
               (l1 + p)^2, (l2 + p)^2};
numeric     = {s -> -3, msq -> 1};
targets     = {j[sunrise, 1, 1, 1, 0, 0]};
precision   = 25;
epsorder    = 3;
