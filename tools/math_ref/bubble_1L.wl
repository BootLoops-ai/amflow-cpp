(* ::Package:: *)

(* Reference computation of the 1-loop bubble using AMFlow + Kira.
   Input mirrors tests/test_layer16_amfsystem.cpp::BubblePipeline_RunsToCompletion :
     loop = {l}, leg = {p}, conservation = {}, replacement = {p^2 -> s},
     propagators = { l^2 - msq, (l - p)^2 - msq },
     numeric = { s -> 5, msq -> 1 },
     target  = j[bubblefam, 1, 1].
   Result is written to tests/data/math_ref/bubble_1L.json in a strict
   machine-parsable form so that C++ tests can ingest it directly. *)

current = If[$FrontEnd === Null, $InputFileName, NotebookFileName[]] // DirectoryName;
amfdir  = FileNameJoin[{current, "..", "..", "reference", "amflow-master"}];
Get[FileNameJoin[{amfdir, "AMFlow.m"}]];

(* Tell AMFlow to use Kira, then override the executable paths set by install.m. *)
SetReductionOptions["IBPReducer" -> "Kira"];
$KiraExecutable   = "/usr/local/bin/kira";
$FermatExecutable = "/usr/share/Ferl7/fer64";

(* Family configuration. *)
AMFlowInfo["Family"]       = bubblefam;
AMFlowInfo["Loop"]         = {l};
AMFlowInfo["Leg"]          = {p};
AMFlowInfo["Conservation"] = {};
AMFlowInfo["Replacement"]  = {p^2 -> s};
AMFlowInfo["Propagator"]   = {l^2 - msq, (l - p)^2 - msq};
AMFlowInfo["Numeric"]      = {s -> 5, msq -> 1};
AMFlowInfo["NThread"]      = 4;

precision = 30;
epsorder  = 4;
target    = {j[bubblefam, 1, 1]};

sol = SolveIntegrals[target, precision, epsorder];
Print["==RAW=="];
Print[sol];

(* Serialize as Laurent expansion in eps. *)
toLaurent[expr_] := Module[{ser, coeffs, leading, n, body, formatted},
  ser = Series[expr /. Rule[__, v_] :> v, {eps, 0, epsorder}];
  If[Head[ser] =!= SeriesData,
    (* Already a number / rational *)
    Return[{{0, ToString[N[expr, precision], InputForm]}}]];
  leading = ser[[4]];
  coeffs  = ser[[3]];
  Table[{leading + n - 1, ToString[N[coeffs[[n]], precision], InputForm]},
        {n, 1, Length[coeffs]}]
];

records = Table[
  Module[{key, val, expr, lst},
    key = ToString[sol[[i, 1]], InputForm];
    expr = sol[[i, 2]];
    lst  = toLaurent[expr];
    {key, lst}], {i, 1, Length[sol]}];

(* Emit a hand-rolled JSON: keep it deterministic, easy to parse. *)
emitNum[s_] := "\"" <> s <> "\"";
emitTerm[t_] := "[" <> ToString[t[[1]]] <> ", " <> emitNum[t[[2]]] <> "]";
emitRecord[r_] := "    {\n      \"integral\": " <> emitNum[r[[1]]] <>
  ",\n      \"laurent\": [" <>
  StringJoin[Riffle[emitTerm /@ r[[2]], ", "]] <> "]\n    }";

json = "{\n  \"family\": \"bubblefam\",\n  \"numeric\": { \"s\": 5, \"msq\": 1 },\n  \"precision\": " <>
  ToString[precision] <> ",\n  \"epsorder\": " <> ToString[epsorder] <>
  ",\n  \"results\": [\n" <>
  StringJoin[Riffle[emitRecord /@ records, ",\n"]] <>
  "\n  ]\n}\n";

outdir = FileNameJoin[{current, "..", "..", "tests", "data", "math_ref"}];
If[!DirectoryQ[outdir], CreateDirectory[outdir]];
outfile = FileNameJoin[{outdir, "bubble_1L.json"}];
WriteString[outfile, json];
Print["==WROTE== ", outfile];

Quit[];
