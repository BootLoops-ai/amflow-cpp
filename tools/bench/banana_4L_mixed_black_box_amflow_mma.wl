(* ::Package:: *)

(*
  4-loop banana with TWO different propagator masses (mAsq, mBsq) in
  alternating A,B,A,B,A pattern across the 5 internal lines.
  Combines the 4-loop topology axis (3.F) with the mixed-mass axis
  (3.I) — neither covered jointly before.
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

AMFlowInfo["Family"] = banana4mix;
AMFlowInfo["Loop"] = {l1, l2, l3, l4};
AMFlowInfo["Leg"] = {p1};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p1^2 -> psq};
AMFlowInfo["Propagator"] = {
  l1^2 - mAsq,
  l2^2 - mBsq,
  l3^2 - mAsq,
  l4^2 - mBsq,
  (l1 + l2 + l3 + l4 - p1)^2 - mAsq,
  (l1 + p1)^2,
  (l2 + p1)^2,
  (l3 + p1)^2,
  (l4 + p1)^2,
  (l1 - l2)^2,
  (l1 - l3)^2,
  (l1 - l4)^2,
  (l2 - l3)^2,
  (l2 - l4)^2
};
AMFlowInfo["Numeric"] = {psq -> -3, mAsq -> 1, mBsq -> 2};
AMFlowInfo["NThread"] = 4;

targets = {j[banana4mix, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "banana_4L_mixed_mma_cache"}]];

Print["==BENCH== mma 4-loop banana mixed-mass BlackBoxAMFlow"];
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
