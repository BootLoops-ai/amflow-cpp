(* ::Package:: *)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
repo = FileNameJoin[{current, "..", ".."}];
Get[FileNameJoin[{repo, "reference", "amflow-master", "AMFlow.m"}]];

SetReductionOptions["IBPReducer" -> "Kira", "DeleteBlackBoxDirectory" -> False];
Kira`$KiraExecutable = "/usr/local/bin/kira";
Kira`$FermatExecutable = "/usr/share/Ferl7/fer64";

SetAMFOptions["ChopPre" -> 20, "RationalizePre" -> 100,
  "LearnXOrder" -> -1, "TestXOrder" -> 5];

AMFlowInfo["Family"] = box1;
AMFlowInfo["Loop"] = {l};
AMFlowInfo["Leg"] = {p1, p2, p3, p4};
AMFlowInfo["Conservation"] = {p4 -> -p1 - p2 - p3};
AMFlowInfo["Replacement"] = {
  p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p4^2 -> 0,
  (p1 + p2)^2 -> s, (p1 + p3)^2 -> t
};
AMFlowInfo["Propagator"] = {
  l^2,
  (l + p1)^2,
  (l + p1 + p2)^2,
  (l + p1 + p2 + p4)^2
};
AMFlowInfo["Numeric"] = {s -> 100, t -> -1};
AMFlowInfo["NThread"] = 4;

targets = {j[box1, 1, 0, 1, 0]};
precision = 10;
epsorder = 3;
cache = Environment["AMF_BENCH_MMA_CACHE"];
If[cache === $Failed || cache === None,
  cache = "/tmp/desolver_bench_box1_mma3"];

Print["==BENCH== mma box1 single SolveIntegrals-equivalent"];
Print["==CACHE== ", cache];

total = AbsoluteTiming[
  {epslist0, workpre, xorder} = GenerateNumericalConfig[precision, epsorder];
  epslist = epslist0 + (4 - $D0)/2;
  Print["==CONFIG== samples=", Length[epslist], " working_pre=", workpre,
    " x_order=", xorder];

  Block[{$WorkingPre = workpre, $XOrder = xorder},
    blackbox = AbsoluteTiming[sol = BlackBoxAMFlow[targets, epslist, cache];][[1]];
    leading = -2*Length[Loop];
    fit = AbsoluteTiming[
      final = Thread[Keys[sol] -> N[
        Normal@Series[FitEps[epslist0, #, leading]& /@ Values[sol],
          {$Eps, 0, epsorder + leading}], precision]];
    ][[1]];
  ];
][[1]];

Print["==BLACKBOX_TIME== ", blackbox];
Print["==FIT_TIME== ", fit];
Print["==TOTAL_TIME== ", total];
Print["==RESULT== ", final // InputForm];

Quit[];
