(* ::Package:: *)

(*
  3-loop banana with 1 massive + 3 massless propagators — Phase 3
  oracle (mixed-mass axis: some massive, some massless).  Same
  banana topology as banana_3loop but with only the first internal
  line carrying msq; the other 3 internal lines are massless.  This
  exercises the scaleless-sub-sector detection and mass-injection
  AMFMode under mixed-mass conditions, neither of which is covered
  by existing benches (banana_3loop is all-equal-mass, tt_2loop_box
  has mostly-massless mixed but only at 2 loops).  Target is the
  corner j[bn3mix, 1,1,1,1, 0,0,0,0,0] at psq = -3, msq = 1,
  eps = 1/1000.
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions[
  "IBPReducer" -> "Kira",
  "BlackBoxRank" -> 0,
  "BlackBoxDot" -> 5,
  "DeleteBlackBoxDirectory" -> False
];
Kira`$KiraExecutable = "/usr/local/bin/kira";
Kira`$FermatExecutable = "/usr/share/Ferl7/fer64";

SetAMFOptions[
  "AMFMode" -> {"Prescription", "Mass", "Propagator"},
  "EndingScheme" -> {"Tradition", "SingleMass"},
  "WorkingPre" -> 120,
  "ChopPre" -> 20,
  "RationalizePre" -> 100,
  "XOrder" -> 200,
  "ExtraXOrder" -> 240,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = bn3mix;
AMFlowInfo["Loop"] = {l1, l2, l3};
AMFlowInfo["Leg"] = {p1};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p1^2 -> psq};
AMFlowInfo["Propagator"] = {
  l1^2 - msq,
  l2^2,
  l3^2,
  (l1 + l2 + l3 - p1)^2,
  (l1 + p1)^2,
  (l2 + p1)^2,
  (l3 + p1)^2,
  (l1 - l2)^2,
  (l1 - l3)^2
};
AMFlowInfo["Numeric"] = {psq -> -3, msq -> 1};
AMFlowInfo["NThread"] = 4;

targets = {
  j[bn3mix, 1, 1, 1, 1, 0, 0, 0, 0, 0]
};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "bn3mix_eps001_mma_cache"}]];

Print["==BENCH== mma 3-loop mixed-mass banana (1m+3 massless) BlackBoxAMFlow sampled"];
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

Quit[];
