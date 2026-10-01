(* ::Package:: *)

(*
  1-loop hexagon — Phase 3 oracle (external-leg axis = 6).
  All-massless internal lines, 6 external on-shell legs (p1..p6 with
  p6 = -(p1+p2+p3+p4+p5) implicit), 6 cyclic propagators.  Extends
  1-loop pentagon to one more leg.  Pair-only Replacement rules
  (lesson from pentagon: no triple-sum identities).  6-pt Mandelstam
  closure: sum_{i<j, i,j in 1..5} s_ij = 0 (since p6^2 = 0).  Chosen
  10 rationals satisfy this:
    s12 -2, s13 +1, s14 +3, s15 +5, s23 -3,
    s24 -2, s25 +1, s34 -5, s35 +4, s45 -2  (sum = 0).
  Target: corner j[hexagon, 1, 1, 1, 1, 1, 1] at eps = 1/1000.
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
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
  "WorkingPre" -> 120,
  "ChopPre" -> 20,
  "RationalizePre" -> 100,
  "XOrder" -> 240,
  "ExtraXOrder" -> 280,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = hexagon;
AMFlowInfo["Loop"] = {l};
AMFlowInfo["Leg"] = {p1, p2, p3, p4, p5};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0, p5^2 -> 0,
  (p1 + p2)^2 -> s12,
  (p1 + p3)^2 -> s13,
  (p1 + p4)^2 -> s14,
  (p1 + p5)^2 -> s15,
  (p2 + p3)^2 -> s23,
  (p2 + p4)^2 -> s24,
  (p2 + p5)^2 -> s25,
  (p3 + p4)^2 -> s34,
  (p3 + p5)^2 -> s35,
  (p4 + p5)^2 -> s45
};
AMFlowInfo["Propagator"] = {
  l^2,
  (l + p1)^2,
  (l + p1 + p2)^2,
  (l + p1 + p2 + p3)^2,
  (l + p1 + p2 + p3 + p4)^2,
  (l + p1 + p2 + p3 + p4 + p5)^2
};
AMFlowInfo["Numeric"] = {
  s12 -> -2, s13 -> 1, s14 -> 3, s15 -> 5,
  s23 -> -3, s24 -> -2, s25 -> 1,
  s34 -> -5, s35 -> 4,
  s45 -> -2
};
AMFlowInfo["NThread"] = 4;

targets = { j[hexagon, 1, 1, 1, 1, 1, 1] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "hexagon_1L_eps001_mma_cache"}]];

Print["==BENCH== mma 1-loop massless hexagon BlackBoxAMFlow sampled"];
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
