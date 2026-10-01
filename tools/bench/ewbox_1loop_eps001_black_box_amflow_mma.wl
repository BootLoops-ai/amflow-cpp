(* ::Package:: *)

(*
  1-loop electroweak box — Phase 3 oracle (multi-invariant axis).
  Family: 4 propagators carrying alternating W/Z masses, 4 massless
  external legs, 4 distinct kinematic invariants {s, t, mWsq, mZsq}.
  Target is the corner j[ewbox, 1, 1, 1, 1] at s = 7, t = -3,
  mWsq = 1, mZsq = 4/3, eps = 1/1000.  Same-MMA-run sampling
  picks up the four 3-propagator subsectors so the oracle covers
  several distinct mass/momentum combinations from one execution.
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
  "WorkingPre" -> 120,
  "ChopPre" -> 20,
  "RationalizePre" -> 100,
  "XOrder" -> 200,
  "ExtraXOrder" -> 240,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = ewbox;
AMFlowInfo["Loop"] = {l};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1 - p2 - p3};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0,
  p2^2 -> 0,
  p3^2 -> 0,
  p4^2 -> 0,
  (p1 + p2)^2 -> s,
  (p1 + p3)^2 -> t
};
AMFlowInfo["Propagator"] = {
  l^2 - mWsq,
  (l + p1)^2 - mZsq,
  (l + p1 + p2)^2 - mWsq,
  (l + p1 + p2 + p3)^2 - mZsq
};
AMFlowInfo["Numeric"] = {s -> 7, t -> -3, mWsq -> 1, mZsq -> 4/3};
AMFlowInfo["NThread"] = 4;

targets = {
  j[ewbox, 1, 1, 1, 1]
};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "ewbox_1loop_eps001_mma_cache"}]];

Print["==BENCH== mma 1-loop electroweak box BlackBoxAMFlow sampled"];
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
    N[Re[val], 40] // InputForm, " :: ", N[Im[val], 40] // InputForm],
  {i, Length[Keys[sol]]}
];

Quit[];
