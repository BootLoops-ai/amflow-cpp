(* ::Package:: *)

(*
  Non-planar 2L crossed-box (xbox) with mass on the CROSS-RUNG
  propagator (l1+l2)^2 - msq.  Currently the suite covers xbox only
  in two massless variants (xbox_2loop_eps001, xbox_2L_eps10).  This
  bench adds a mass on the specific propagator that joins both loops
  in the non-planar (crossed) configuration — the analogue of
  doublebox_diagmass for the non-planar topology.

  Combines NON-PLANAR + MASS code paths.  s = -3, t = -1 (Euclidean),
  msq = 1.  eps = 1/1000.
  Target: corner j[xboxCrossmass, 1,1,1,1,1,1,1,0,0].
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions[
  "IBPReducer" -> "Kira",
  "BlackBoxRank" -> 0,
  "BlackBoxDot" -> 1,
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

AMFlowInfo["Family"] = xboxCrossmass;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1 - p2 - p3};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
  (p1 + p2)^2 -> s,
  (p1 + p3)^2 -> t
};
AMFlowInfo["Propagator"] = {
  l1^2,
  (l1 + p1)^2,
  (l1 + p1 + p2)^2,
  l2^2,
  (l2 - p3)^2,
  (l2 + p1 + p2)^2,
  (l1 + l2)^2 - msq,
  (l1 - p3)^2,
  (l2 + p1)^2
};
AMFlowInfo["Numeric"] = {
  s -> -3, t -> -1, msq -> 1
};
AMFlowInfo["NThread"] = 4;

targets = { j[xboxCrossmass, 1, 1, 1, 1, 1, 1, 1, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "xbox_crossmass_2L_eps001_mma_cache"}]];

Print["==BENCH== mma 2-loop non-planar xbox with mass on cross-rung"];
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
