(* ::Package:: *)

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
  "XOrder" -> 200,
  "ExtraXOrder" -> 240,
  "LearnXOrder" -> -1,
  "TestXOrder" -> 5
];

AMFlowInfo["Family"] = bn34m;
AMFlowInfo["Loop"] = {l1, l2, l3};
AMFlowInfo["Leg"] = {p1};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"] = {p1^2 -> psq};
AMFlowInfo["Propagator"] = {
  l1^2 - mAsq,
  l2^2 - mBsq,
  l3^2 - mCsq,
  (l1 + l2 + l3 - p1)^2 - mDsq,
  (l1 + p1)^2,
  (l2 + p1)^2,
  (l3 + p1)^2,
  (l1 - l2)^2,
  (l1 - l3)^2
};
AMFlowInfo["Numeric"] = {psq -> -3, mAsq -> 1, mBsq -> 4, mCsq -> 9, mDsq -> 16};
AMFlowInfo["NThread"] = 4;

targets = {
  j[bn34m, 1, 1, 1, 1, 0, 0, 0, 0, 0]
};
epslist = {1/1000};
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = FileNameJoin[{repo, "tools", "bench", "mma_refs",
    "bn3_4mass_3L_eps001_mma_cache"}]];

Print["==BENCH== mma bn3_4mass 3L 4-distinct-mass BlackBoxAMFlow sampled"];
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

(* Dump the AMFSystem-1 master-integral solution table for sampled
   parity (`sampled_values_from_same_mma_run` in the reference JSON).
   `solution` is MMA's `Put`-format file of {master -> {eps_values...}}. *)
sysSol = Get[FileNameJoin[{cache, "1", "solution"}]];
If[Head[sysSol] === List,
  Do[
    key = sysSol[[i, 1]];
    val = sysSol[[i, 2, 1]];
    Print["==SAMPLED== ", key // InputForm, " :: ",
      N[Re[val], 50] // InputForm, " :: ", N[Im[val], 50] // InputForm],
    {i, Length[sysSol]}
  ]];

Quit[];
