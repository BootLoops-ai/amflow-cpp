(* ::Package:: *)

(*
  Doublebox 2L 4-leg with mass on the INTER-LOOP propagator
  (l1 - l2)^2 - msq.  Novel mass placement axis: existing 2L mass
  benches (bn3mix, wzbox, doublebox2m, doublebox_sv) put mass either
  within one loop or as a "block" pattern; none places mass on the
  rung propagator shared between the two loops.  This tests the
  η-injection + AMFMode dispatch on the inter-loop propagator code
  path.

  Same Mandelstam kinematics as doublebox_sv:
    s = -3, t = -1 (Euclidean), msq = 1.
  Target: corner j[doubleboxDiagmass, 1,1,1,1,1,1,1,0,0] at eps = 1/1000.
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

AMFlowInfo["Family"] = doubleboxDiagmass;
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
  (l2 + p1 + p2)^2,
  (l2 + p1 + p2 + p4)^2,
  l2^2,
  (l1 - l2)^2 - msq,
  (l1 + p3)^2,
  (l2 + p1)^2
};
AMFlowInfo["Numeric"] = {
  s -> -3, t -> -1, msq -> 1
};
AMFlowInfo["NThread"] = 4;

targets = { j[doubleboxDiagmass, 1, 1, 1, 1, 1, 1, 1, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "doubleboxDiagmass_2L_eps001_mma_cache"}]];

Print["==BENCH== mma 2-loop doublebox with diagonal-rung mass BlackBoxAMFlow"];
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
