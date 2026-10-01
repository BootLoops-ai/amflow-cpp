(* ::Package:: *)

(*
  2L 3-leg vertex (vtx2 topology) with THREE distinct masses placed
  on three different propagators (mAsq on l1^2, mBsq on
  (l1+p1+p2)^2, mCsq on l2^2).  Different from the all-massless
  vtx2_2loop_vertex_eps001 oracle: this variant exercises the
  multi-mass code path at 2L 3-leg.

  Hierarchical masses (1, 4, 9) discriminate the SingleMass
  scheme dispatch.  s = -1 (Euclidean off-shell).  eps = 1/1000.
  Target: corner j[vtx23m, 1, 1, 1, 1, 1, 1, 0].
*)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions[
  "IBPReducer" -> "Kira",
  "BlackBoxRank" -> 2,
  "BlackBoxDot" -> 0,
  "DeleteBlackBoxDirectory" -> False
];
Kira`$KiraExecutable = "/usr/local/bin/kira";
Kira`$FermatExecutable = "/usr/share/Ferl7/fer64";

SetAMFOptions[
  "AMFMode" -> {"Prescription", "Mass", "Propagator"},
  "EndingScheme" -> {"Tradition", "SingleMass"},
  "WorkingPre" -> 120, "ChopPre" -> 20, "RationalizePre" -> 100,
  "XOrder" -> 200, "ExtraXOrder" -> 240,
  "LearnXOrder" -> -1, "TestXOrder" -> 5
];

AMFlowInfo["Family"] = vtx23m;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0, p2^2 -> 0, (p1 + p2)^2 -> s
};
AMFlowInfo["Propagator"] = {
  l1^2 - mAsq,
  (l1 + p1)^2,
  (l1 + p1 + p2)^2 - mBsq,
  l2^2 - mCsq,
  (l2 + p1 + p2)^2,
  (l1 - l2)^2,
  (l2 + p1)^2
};
AMFlowInfo["Numeric"] = {
  s -> -1, mAsq -> 1, mBsq -> 4, mCsq -> 9
};
AMFlowInfo["NThread"] = 4;

targets = { j[vtx23m, 1, 1, 1, 1, 1, 1, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "vtx2_2L_3mass_eps001_mma_cache"}]];

Print["==BENCH== mma 2L vertex with 3 distinct masses"];
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
