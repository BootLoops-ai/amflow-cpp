(* ::Package:: *)

(*
  η_c → γγ 3-loop Mercedes-Benz form factor — STAGED MASS PROBE (mass1)
  Rim-insertion topology (3-fold symmetric, all 3-valent / φ^3).

  In this variant, ONLY propagator 1 (l1² - msq) is massive (placeholder
  c-quark). All other 5 rim segments are massless along with the spokes.
  Purpose: smallest deviation from the validated massless oracle, used to
  measure (n_massive, dot_needed, runtime) curve for pseudo-master onset.

  Kinematics (Euclidean):
    p1 + p2 + p3 = 0 (all-incoming), p1^2 = ssq = -3, p2^2 = p3^2 = 0.
    msq = 1.
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions[
  "IBPReducer" -> "Kira",
  "BlackBoxRank" -> 0,
  "BlackBoxDot" -> 3,
  "DeleteBlackBoxDirectory" -> False
];
Kira`$KiraExecutable = "/usr/local/bin/kira";
Kira`$FermatExecutable = "/usr/share/Ferl7/fer64";

SetAMFOptions[
  "AMFMode" -> {"Prescription", "Mass", "Propagator"},
  "EndingScheme" -> {"Tradition", "SingleMass"},
  "WorkingPre" -> 120, "ChopPre" -> 20, "RationalizePre" -> 100,
  "XOrder" -> 240, "ExtraXOrder" -> 280,
  "LearnXOrder" -> -1, "TestXOrder" -> 5
];

AMFlowInfo["Family"] = mercedesM1;
AMFlowInfo["Loop"] = {l1, l2, l3};
AMFlowInfo["Leg"] = {p1, p2, p3};
AMFlowInfo["Conservation"] = {p3 -> -p1 - p2};
AMFlowInfo["Replacement"] = {
  p1^2 -> ssq,
  p2^2 -> 0,
  (p1 + p2)^2 -> 0
};
AMFlowInfo["Propagator"] = {
  l1^2 - msq,
  (l1 + p1)^2,
  l2^2,
  (l2 + p2)^2,
  l3^2,
  (l3 + p3)^2,
  (l2 - l1 - p1)^2,
  (l3 - l2 - p2)^2,
  (l1 - l3 - p3)^2,
  (l1 + p2)^2,
  (l2 + p3)^2,
  (l3 + p1)^2
};
AMFlowInfo["Numeric"] = {
  ssq -> -3,
  msq -> 1
};
AMFlowInfo["NThread"] = 4;

targets = { j[mercedesM1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "eta_c_form_factor", "mma_refs",
    "eta_c_3L_Mercedes_mass1_euclid_eps001_mma_cache"}]];

Print["==BENCH== mma eta_c 3L Mercedes mass1 staged (1 c-quark, prop 1)."];
Print["==CACHE== ", cache];
Print["==EPS== ", epslist // InputForm];

elapsed = AbsoluteTiming[
  sol = BlackBoxAMFlow[targets, epslist, cache];
][[1]];

Print["==TIME== ", elapsed];
Print["==RESULT_INPUTFORM== ", sol // InputForm];

Do[
  key = Keys[sol][[i]];
  val = Values[sol][[i, 1]];
  Print["==VALUE== ", key // InputForm, " :: ",
    N[Re[val], 50] // InputForm, " :: ", N[Im[val], 50] // InputForm],
  {i, Length[Keys[sol]]}
];

sysSol = Get[FileNameJoin[{cache, "1", "solution"}]];
If[Head[sysSol] === List,
  Do[
    key = sysSol[[i, 1]];
    val = sysSol[[i, 2, 1]];
    Print["==SAMPLED== ", key // InputForm, " :: ",
      N[Re[val], 50] // InputForm, " :: ", N[Im[val], 50] // InputForm],
    {i, Length[sysSol]}
  ]];

Quit[];
