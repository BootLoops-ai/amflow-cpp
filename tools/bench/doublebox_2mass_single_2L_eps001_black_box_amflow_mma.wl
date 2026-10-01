(* ::Package:: *)

(*
  Doublebox 2L 4-leg with SPARSE 2-mass placement: mAsq on l1^2
  (prop 0 only) and mBsq on l2^2 (prop 5 only).  Only 2 of the 9
  propagators carry mass.  Different from:
    - doublebox2m_eps001 (D12): mAsq on l1's prop 0 AND l2's prop 3
      (cross-loop SAME mass interleaved); mBsq on l1's prop 1 AND
      l2's prop 5.
    - doublebox_blockmass_2L: mAsq on l1's two rails (props 0+1) and
      mBsq on l2's two rails (props 3+5); BLOCK pattern, 4 massive
      props.
    - doublebox_diagmass_2L: msq on (l1-l2) inter-loop propagator.
    - doublebox_sv: all-massless.

  This sparse-2mass pattern (one mass per loop on the corner
  propagator only) is the simplest non-trivial 2-mass variation and
  exercises a different region-enumeration path than block or
  interleaved.

  s = -3, t = -1 (Euclidean), mAsq = 1, mBsq = 4.  eps = 1/1000.
  Target: corner j[doubleboxSingle2m, 1,1,1,1,1,1,1,0,0].
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

AMFlowInfo["Family"] = doubleboxSingle2m;
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
  (l1 + p1)^2,
  (l1 + p1 + p2)^2,
  (l2 + p1 + p2)^2,
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

targets = { j[doubleboxSingle2m, 1, 1, 1, 1, 1, 1, 1, 0, 0] };
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "doublebox_2mass_single_2L_eps001_mma_cache"}]];

Print["==BENCH== mma 2L doublebox with sparse 2-mass (one per loop)"];
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
