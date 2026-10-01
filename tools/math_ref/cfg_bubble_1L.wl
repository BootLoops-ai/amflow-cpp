(* 1-loop massive bubble matching tests/test_layer16_amfsystem.cpp.
   We push precision/epsorder high enough so that evaluating the
   truncated Laurent at eps=1/100 is accurate to <1e-12. *)
name        = "bubble_1L";
family      = bubblefam;
loops       = {l};
legs        = {p};
conservation = {};
replacement = {p^2 -> s};
propagators = {l^2 - msq, (l - p)^2 - msq};
numeric     = {s -> 5, msq -> 1};
targets     = {j[bubblefam, 1, 1], j[bubblefam, 1, 0], j[bubblefam, 0, 1]};
precision   = 50;
epsorder    = 10;
