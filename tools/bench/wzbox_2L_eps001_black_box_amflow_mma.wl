(* ::Package:: *)

(*
  2-loop W/Z box — Phase 3 oracle (combined-axis: L=2 + 6 invariants +
  alternating heavy internal masses + cross-threshold).  4 external
  massless legs (s, t kinematics), 2 loops, 9 propagators with
  alternating mW²/mZ² masses on the 6 main rails and 3 massless
  ISP/coupling propagators.  Extends the 1-loop ewbox (W/Z alternating)
  to 2 loops.  Numerics: s = 7 (> 4 mW² = 4, above W-pair threshold ⇒
  imaginary parts), t = -3, mWsq = 1, mZsq = 4/3, eps = 1/1000.
  BlackBoxDot = 5 to avoid the upstream "inconsistent masters" abort
  pattern that plagued 5L banana / 3L double-box at low dot.
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
  "XOrder" -> 240,
  "ExtraXOrder" -> 280,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = wzbox2L;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1 - p2 - p3};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
  (p1 + p2)^2 -> s,
  (p1 + p3)^2 -> t
};
AMFlowInfo["Propagator"] = {
  l1^2 - mWsq,
  (l1 + p1)^2 - mZsq,
  (l1 + p1 + p2)^2 - mWsq,
  l2^2 - mZsq,
  (l2 - p3)^2 - mWsq,
  (l2 + p1 + p2)^2 - mZsq,
  (l1 + l2)^2,
  (l1 - p3)^2,
  (l2 + p1)^2
};
AMFlowInfo["Numeric"] = {s -> 7, t -> -3, mWsq -> 1, mZsq -> 4/3};
AMFlowInfo["NThread"] = 4;

targets = { j[wzbox2L, 1, 1, 1, 1, 1, 1, 1, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "wzbox_2L_eps001_mma_cache"}]];

Print["==BENCH== mma 2-loop W/Z box BlackBoxAMFlow sampled"];
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
