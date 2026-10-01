(* ::Package:: *)

(*
  Doublebox 2L 4-leg with BLOCK 2-mass pattern (mAsq on l1's two
  rails — propagators 0 and 1; mBsq on l2's two rails — propagators
  3 and 5).  Direct topological complement to doublebox2m (D12
  regression guard, which has the same family with INTERLEAVED
  cross-loop mass placement).

  The audit's D12 narrative argued that the block placement gives a
  top-sector diffeq matrix with real-axis-only poles (denominators
  like (eta+1)), whereas the interleaved placement produces complex
  pole pair η^2 + η + 12 that tightens the Frobenius convergence
  radius.  This bench locks the block case at full default precision
  and gives a side-by-side comparison standard for the interleaved
  precision-sensitive case.

  Same kinematic numerics as doublebox_sv / doublebox2m: s=-3, t=-1,
  mAsq=1, mBsq=4 (Euclidean).  eps = 1/1000.
  Target: corner j[doubleboxBlockmass, 1,1,1,1,1,1,1,0,0].
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

AMFlowInfo["Family"] = doubleboxBlockmass;
AMFlowInfo["Loop"] = {l1, l2};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1 - p2 - p3};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
  (p1 + p2)^2 -> s,
  (p1 + p3)^2 -> t
};
AMFlowInfo["Propagator"] = {
  l1^2 - mAsq,
  (l1 + p1)^2 - mAsq,
  (l1 + p1 + p2)^2,
  (l2 + p1 + p2)^2 - mBsq,
  (l2 + p1 + p2 + p4)^2,
  l2^2 - mBsq,
  (l1 - l2)^2,
  (l1 + p3)^2,
  (l2 + p1)^2
};
AMFlowInfo["Numeric"] = {
  s -> -3, t -> -1, mAsq -> 1, mBsq -> 4
};
AMFlowInfo["NThread"] = 4;

targets = { j[doubleboxBlockmass, 1, 1, 1, 1, 1, 1, 1, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "doublebox_blockmass_2L_eps001_mma_cache"}]];

Print["==BENCH== mma 2-loop doublebox with block 2-mass scheme"];
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
