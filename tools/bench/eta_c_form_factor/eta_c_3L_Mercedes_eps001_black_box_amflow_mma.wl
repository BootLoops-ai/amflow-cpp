(* ::Package:: *)

(*
  η_c → γγ style 3-loop Mercedes-Benz form factor
  Rim-insertion topology (3-fold symmetric, all 3-valent / φ^3).

  Geometry:
    - Central vertex O (3 spokes to A, B, C)
    - Outer rim: A-D-B-E-C-F-A (6 segments, externals at D, E, F)
    - 9 internal propagators (3 spokes + 6 rim segments), 7 vertices, 3 loops.
    - External legs: p1 at D, p2 at E, p3 at F.

  Kinematics: 1 -> 2 decay
    - p1 + p2 + p3 = 0 (all-incoming)
    - p1^2 = ssq (incoming massive, off-shell)
    - p2^2 = p3^2 = 0 (outgoing massless)
    - Single scale: ssq = 1.
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

AMFlowInfo["Family"] = mercedes;
AMFlowInfo["Loop"] = {l1, l2, l3};
AMFlowInfo["Leg"] = {p1, p2, p3};
AMFlowInfo["Conservation"] = {p3 -> -p1 - p2};
AMFlowInfo["Replacement"] = {
  p1^2 -> ssq,
  p2^2 -> 0,
  (p1 + p2)^2 -> 0
};
AMFlowInfo["Propagator"] = {
  l1^2,
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
  ssq -> 1
};
AMFlowInfo["NThread"] = 4;

targets = { j[mercedes, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "eta_c_form_factor", "mma_refs",
    "eta_c_3L_Mercedes_eps001_mma_cache"}]];

Print["==BENCH== mma eta_c 3L Mercedes-Benz form factor (rim-insertion, 9 prop + 3 ISP)."];
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
