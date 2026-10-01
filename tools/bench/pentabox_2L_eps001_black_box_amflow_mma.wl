(* ::Package:: *)

(*
  2-loop pentabox — Phase 3 oracle.  Combines L=2 (existing) with 5
  external legs (pentagon already at 1L) — first 2L 5-point bench.
  All-massless internals, 4 reduced legs after conservation
  p5 = -(p1+p2+p3+p4), 8 main propagators + 3 ISPs forming the 2-loop
  SP basis of size 11.  Numerics: same 5-pt Mandelstam basis as 1L
  pentagon (s12=-2, s23=-3, s34=-5, s45=-7, s51=-11; derived s13=-2,
  s14=+15, s24=-3, s25=+8, s35=+10), eps = 1/1000.
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
  "WorkingPre" -> 120, "ChopPre" -> 20, "RationalizePre" -> 100,
  "XOrder" -> 240, "ExtraXOrder" -> 280,
  "LearnXOrder" -> -1, "TestXOrder" -> 5
];

AMFlowInfo["Family"] = pentabox2L;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
  (p1 + p2)^2 -> s12,
  (p1 + p3)^2 -> s13,
  (p1 + p4)^2 -> s14,
  (p2 + p3)^2 -> s23,
  (p2 + p4)^2 -> s24,
  (p3 + p4)^2 -> s34
};
AMFlowInfo["Propagator"] = {
  l1^2,
  (l1 + p1)^2,
  (l1 + p1 + p2)^2,
  (l1 + p1 + p2 + p3)^2,
  l2^2,
  (l2 - p4)^2,
  (l1 + l2)^2,
  (l2 + p1 + p2)^2,
  (l1 - p4)^2,
  (l2 + p1)^2,
  (l2 - p3)^2
};
AMFlowInfo["Numeric"] = {
  s12 -> -2, s23 -> -3, s34 -> -5,
  s13 -> -2, s14 -> 15, s24 -> -3
};
AMFlowInfo["NThread"] = 4;

targets = { j[pentabox2L, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "pentabox_2L_eps001_mma_cache"}]];

Print["==BENCH== mma 2-loop pentabox BlackBoxAMFlow sampled"];
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
