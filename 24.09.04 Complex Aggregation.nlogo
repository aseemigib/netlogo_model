turtles-own [
  folded? ; TRUE at setup; FALSE after unfolding
  leader  ; used to coordinate movements of aggregates
  species
  weight
  deltaG Keq
  k
]

globals [
  h
  scaling_factor
  scribe0 scribe1 ; for the monitors
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;BEGINNING OF SETUP COMMANDS;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all no-display
  reset-timer

  if invalid_inputs? [stop]

  setup_globals
  setup_proteins

  update_locals
  update_aesthetics

  reset-ticks display
end

to-report invalid_inputs?
  if sum props != 1 [print "error: proportions must sum to 1" report TRUE]
  if length (runresult weights) != length props [print "error: inputs must be of same length" report TRUE]
  report FALSE
end

to setup_globals ; these parameters differ from protein to protein and could be added in an extension
  set h 25
  set scaling_factor .004
end

to setup_proteins
  create-turtles population [
    set leader self
    set shape "circle"
    set folded? TRUE
    setxy random-xcor random-ycor
    set species who_to_species who
  ]
  set scribe0 one-of turtles with [species = 0]
  set scribe1 one-of turtles with [species = 1]
end

to update_locals
  ask turtles [
    set weight item species runresult weights ;ifelse-value species = 0 [weight0] [weight1]
    set size weight / 20

    set k ( 0.02 * (weight / .11))
    set deltaG (scaling_factor * ((temperature - h) ^ 2) - k)
    set Keq (e ^ ((- deltaG * 1000) / (8.314 * (temperature + 273.15)) ))
  ]
end

to update_aesthetics
  if behaviorspace-run-number = 0 [ ; if not in behaviorspace
    ask turtles [
      ifelse folded?
      [ set color [0 255 0 100] ]
      [ set color [255 0 0 100] ]

      if any? in-link-neighbors [
        set color [0 0 255 100]
        if not soluble? self [set color [255 255 255 100]]
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;END OF SETUP COMMANDS;;;;;;;;;;;;
;;;;;;;;;;;;BEGINNING OF DYNAMICS;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go
  move_and_fold

  agg_phase
  disagg_phase

  update_aesthetics
  tick
end

;;;;;;;;;;;;;;OBSERVER COMMANDS;;;;;;;;;;;;;;

to move_and_fold ; called by GO
  ask turtles with [ leader = self ] [
    move

    if not any? in-link-neighbors [ ; if disaggregated
      ifelse folded?
      [ unfold-forward ]
      [ fold-backward ]
    ]
  ]
end

to agg_phase ; called by GO
  ask turtles with [ leader = self ] [
    if not folded? [ aggr ]
  ]
end

to disagg_phase ; called by GO
  ask turtles [
    if probability? disagg_rate [
      let neighborhood link-neighbors

      ask my-links [ die ]
      set leader self

      if any? neighborhood [
        ask neighborhood [
          set leader self
          ask link-neighbors [ merge ]
        ]
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;AGENT  COMMANDS;;;;;;;;;;;;;;;

to move ; turtle procedure
  let agg aggregate self ; precompute this now
  let wt sum ([weight] of agg)

  let head random 360
  let dist (temperature + 273.15) / wt
  ask agg [
    set heading head
    fd dist
  ]
end

to unfold-forward  ; turtle procedure
  let normalizer max (list (Keq) (1 / Keq))

  let chance (1 / Keq) / normalizer
  if probability? (.1 * chance) [set folded? FALSE]
end

to fold-backward  ; turtle procedure
  let normalizer max (list (Keq) (1 / Keq))

  let chance (Keq) / normalizer
  if probability? (.1 * chance) [set folded? TRUE]
end

to aggr ; turtle procedure
  let candidates other turtles in-radius stickiness with [not folded?]
  foreach [self] of candidates [ x ->
    if probability? agg_rate [
      create-link-with x
      ask x [ merge ]
    ]
  ]
end

to merge ; turtle procedure
         ; This is a recursive process

  ; First this protein merges to whomever it was called by.
  set leader [ leader ] of myself

  ; Then neighboring nodes merge to this protein.
  ask link-neighbors with [ leader != [leader] of myself ] [ merge ]
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;END OF MOTOR COMMANDS;;;;;;;;;;;;
;;;;;;BEGINNING OF REPORTER DEFINITIONS;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to-report unf_conc
  report (count turtles with [not folded?]) / (count turtles)
end

to-report fold_conc
  report (count turtles with [folded?]) / (count turtles)
end

to-report total_conc
  report (count turtles)
end

to-report fold0
  let prot0 turtles with [species = 0]
  report count prot0 with [folded?] / count prot0
end

to-report fold1
  let prot0 turtles with [species = 1]
  report count prot0 with [folded?] / count prot0
end

to-report aggregate [turt]
  let ldr [leader] of turt
  report turtles with [leader = ldr]
end

to-report scale [turt]
  report count aggregate turt
end

to-report makeup [turt]
  let agg aggregate turt
  let prot0 agg with [species = 0]
  ifelse any? prot0 [report count prot0 / count agg] [report 0]
end
to-report soluble? [turt]
  ; this reporter is coded in a clunky manner
  ; this is due to the cost of computing the scale
  ; which rises linearly with the population
  ; this reporter is called on every turtle on every tick
  ; which would then rise superlinearly with population

  ; there is also an implicit parameter here
  ; which reflects the arbitrary cut-off of 2

  let neighborhood [link-neighbors] of turt

  if not any? neighborhood [ report true ] ; if no neighbors

  ifelse count neighborhood > 1 [
    report false
  ] [ ; else if multiple neighbors
    report (count [link-neighbors] of one-of neighborhood > 1)
  ]
end

to-report soluble_conc
  report count turtles with [ soluble? self ] / count turtles
end

to-report freq_dist [dist]
  let sorted_dist reverse sort dist
  let values remove-duplicates sorted_dist

  let result []

  foreach values [ value ->
    let freq length (filter [x -> x = value] sorted_dist)
    let new_entree list value freq
    set result insert-item 0 result new_entree
  ]

  report result
end

to-report probability? [ p ] ; reports TRUE with probability p; else reports FALSE
  ifelse (random-float 1 < p)
  [ report TRUE ]
  [ report FALSE ]
end

to-report props
  report runresult proportions
end

to-report who_to_species [w]
  let prop 0
  foreach (range length props) [s ->
    set prop prop + (item s props)
    if (w / population) < prop [report s]
  ]
end

to-report ks
  let l []
  foreach (range length props) [ s ->
    set l lput precision ([k] of one-of turtles with [species = s]) 2 l
  ]
  report l
end

to-report deltaGs
  let l []
  foreach (range length props) [ s ->
    set l lput precision ([deltaG] of one-of turtles with [species = s]) 2 l
  ]
  report l
end

to-report Keqs
  let l []
  foreach (range length props) [ s ->
    set l lput precision ([Keq] of one-of turtles with [species = s]) 2 l
  ]
  report l
end

to-report folded_concs
  let l []
  foreach (range length props) [ s ->
    let candidates turtles with [species = s]
    let targets candidates with [folded?]
    let fc ifelse-value any? targets [count targets / count candidates] [0]
    set l lput precision fc 2 l
  ]
  report l
end

to-report soluble_concs
  let l []
  foreach (range length props) [ s ->
    let candidates turtles with [species = s]
    let targets candidates with [soluble? self]
    let fc ifelse-value any? targets [count targets / count candidates] [0]
    set l lput precision fc 2 l
  ]
  report l
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;END OF COMMANDS;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@#$#@#$#@
GRAPHICS-WINDOW
689
17
1290
619
-1
-1
2.95025
1
12
1
1
1
0
1
1
1
-100
100
-100
100
1
1
1
ticks
30.0

BUTTON
14
46
77
79
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
79
46
142
79
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
460
56
590
89
Temperature
Temperature
0
100
80.0
1
1
C
HORIZONTAL

SLIDER
14
10
206
43
Population
Population
0
10000
2000.0
100
1
NIL
HORIZONTAL

PLOT
9
177
189
337
Folding Dynamics
ticks
Percent of Protein
0.0
10.0
0.0
100.0
true
true
"" ""
PENS
"Unfolded" 1.0 0 -2674135 true "" "plot unf_conc * 100"
"Folded" 1.0 0 -13840069 true "" "plot fold_conc * 100"

BUTTON
145
46
211
79
go once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
405
422
483
467
Soluble %
soluble_conc
3
1
11

PLOT
192
177
371
336
% Soluble over Time
ticks
Percent of Protein
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"Soluble Conc" 1.0 0 -16777216 true "" "plot soluble_conc"

PLOT
35
343
348
463
Aggregate Scales by Frequency
Weight of Aggregate
Frequency
2.0
50.0
0.0
20.0
true
false
"" ""
PENS
"default" 3.0 1 -16777216 true "" "histogram [scale self] of turtles with [leader = self and not folded?]"

SLIDER
123
100
238
133
disagg_rate
disagg_rate
0
0.05
0.02
.005
1
NIL
HORIZONTAL

BUTTON
10
135
154
168
turn off aggregation
set agg_rate 0\nset disagg_rate 0\n;ask links [die]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
160
136
271
169
kill aggregates
ask links [die]\nask turtles [set leader self]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
8
100
116
133
agg_rate
agg_rate
0
1
0.75
.01
1
NIL
HORIZONTAL

SLIDER
243
101
370
134
Stickiness
Stickiness
0
2
1.3
0.1
1
NIL
HORIZONTAL

MONITOR
405
373
482
418
Folded %
fold_conc
3
1
11

PLOT
12
475
331
675
Species-wise aggregates
Aggregate Size
Frequency of Occurance
2.0
60.0
0.0
20.0
true
true
"" ""
PENS
"species 0" 1.0 1 -955883 true "" "histogram [scale self] of turtles with [species = 0 and not folded?]"
"species 1" 1.0 1 -13345367 true "" "histogram [scale self] of turtles with [species = 1 and not folded?]"
"speices 2" 1.0 1 -14439633 true "" "histogram [scale self] of turtles with [species = 2 and not folded?]"

PLOT
334
476
637
674
Aggregates - Composition by Size
Aggregate Size
Percent Species 0
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 2 -16777216 true "" "clear-plot\nforeach sort turtles with [leader = self and not folded?] [agg ->\nlet x (scale agg) + random-normal 0 .1\nlet y (makeup agg) + random-normal 0 .01\nplotxy x y\n]"

INPUTBOX
389
106
635
166
Weights
[60 30]
1
0
String

INPUTBOX
389
174
667
234
Proportions
[0.5 0.5]
1
0
String

MONITOR
523
241
632
286
NIL
deltaGs
17
1
11

MONITOR
406
241
516
286
NIL
ks
17
1
11

MONITOR
408
291
515
336
NIL
Keqs
17
1
11

MONITOR
486
373
591
418
Folded %s
folded_concs
17
1
11

MONITOR
486
422
591
467
Soluble %s
soluble_concs
17
1
11

BUTTON
684
640
790
673
time 100 ticks
let x timer\nrepeat 100 [go]\nlet y 100 / (timer - x)\nprint word y \" ticks per second\"
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## WHAT IS IT?

This is an oversimplified 2D kinetic model of protein thermal aggregation. It is created on the basis of experimental results obtained in thermal aggregation experiments carried out for single protein species. The model extrapolates this single-protein thermal aggregation behaviour to multiple protein species aggregating together as a mixture.

The model can be used to-

1. Get a general feel of protein thermal aggregation
2. Check how thermal aggregation behaviour of protein changes by changing values of different parameters in the code.
3. Compare simulation results with experimental results of protein thermal aggregation using various protein systems
4. Generate new hypotheses about thermal aggregation of complex protein mixtures

Assumptions used in the model (for parsimony)-

1. Proteins exist in one of two states, either folded or unfolded
2. Only unfolded states can bonded to other proteins and form aggregations, and proteins bonded to other proteins cannot fold
3. Thermal aggregation is a reversible process, as proteins can leave the aggregate
4. All proteins within a species will have equal molecular weight
5. The extent of thermal folding-unfolding of protein is a function of its molecular weight
6. All protein species under inspection will form and lose bonds with the same rates
7. Aggregates will be mixed (many different species can be a part of one single aggregate)
8. Proteins can only form links with physically proximate aggregates.

## THEORY BEHIND THE MODEL

![Description of the model](https://lh6.googleusercontent.com/n2Nvr8Ta-uZbwNCGRu8BaiM5OWa5RM8PtxfVnk7A65F55tqHYrun9XDpqk1E-QiEf5jWUOT6FxF6xVAhOhH59L6dNa6vZOmBKs1Dkdo5KfOTAQYZtc4M9wjl4fAxHw4XyA)

Folded to unfolded state transition is driven by elevated temperature, and unfolded to aggregate transition is driven by intermolecular collisions.

Thermodynamics of protein folding-

According to Ganesh, C., et al. "Prediction of the maximal stability temperature of
monomeric globular proteins solely from amino acid sequence." FEBS letters 454.1-2 (1999): 31-36; Gibbs free energy of protein folding (delta G) is directly proportional to protein molecular weight.

Keq = e^(-delta G/RT)

Keq: protein folding equilibrium constant, e: Euler's number, R: universal gas constant, T: temperature

Keq = k1/k-1

This means that the proportion of folded:unfolded protein depends on protein molecular weight and temperature, which are user-defined.

Based on user-defined k_agg, k_disagg and Stickiness, unfolded protein aggregates.

## HOW TO USE IT

Setting up the model:

Choosing desired global and individual agent parameters is the first step. Population, Temperature, Aggregation rate (agg_rate), Disaggregation rate (disagg_rate) and Stickiness are global parameters. On the other hand, Molecular Weights (Weights) of individual protein species and their Relative Proportions (Proportions) are individual agent variables.

First we start the model with folded proteins. When we press 'GO', the folded proteins (green) will start unfolding (red). The rate and extent for unfolding and folding for each individual protein species depends on temperature and their respective molecular weights.

The rates of folding and unfolding naturally emerge from an equilibrium constant (Keq)
which depends on Gibbs free energy of folding, which in turn depends on temperature and molecular weights of protein species.

The output is in the form of two lines in the accompanying plot- native protein (green) and unfolded protein (red).

The model reasonably simulates equilibrium behaviour of protein unfolding process at a
particular temperature. The model is thermodynamically accurate, as it gives value of Keq as 1 when the value of Gibbs free energy change is 0.

For more details, please go to https://sites.google.com/view/simulation-aseem/home
and download 'supplementary information'.

## RELATED MODELS 

This code for this model is based on the Simple Kinetics 2 model from the Netlogo Models Library, in Chemistry and Physics.

Stieff, M. and Wilensky, U. (2001). NetLogo Simple Kinetics 2 model.
http://ccl.northwestern.edu/netlogo/models/SimpleKinetics2.
Center for Connected Learning and Computer-Based Modeling, 
Northwestern University, Evanston, IL.

## CREDITS

This model is a product of a collaboration between Aseem Chaphalkar, Leo Niehorster-Cook, and Dr. Kausik Chakraborty. This project was funded by CSIR Empower grant (MLP2103) and DBT-Wellcome Trust India-Alliance grant IA/S/21/1/505587, both awarded to Dr. Kausik Chakraborty.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
setup
repeat 450 [ go ]
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="individual_conc" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="300"/>
    <metric>soluble_conc</metric>
    <enumeratedValueSet variable="proportion1">
      <value value="0"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="weight0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="weight1">
      <value value="64"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="linkwise?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temperature">
      <value value="25"/>
      <value value="45"/>
      <value value="55"/>
      <value value="65"/>
      <value value="75"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calculate_dynamically?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stickiness">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population">
      <value value="3000"/>
      <value value="4000"/>
      <value value="6000"/>
      <value value="8000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disagg_rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agg_rate">
      <value value="0.75"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="diff_proportion" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="300"/>
    <metric>soluble_conc</metric>
    <enumeratedValueSet variable="proportion1">
      <value value="0.2"/>
      <value value="0.5"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="weight0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="weight1">
      <value value="64"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="linkwise?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temperature">
      <value value="25"/>
      <value value="45"/>
      <value value="55"/>
      <value value="65"/>
      <value value="75"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calculate_dynamically?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stickiness">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population">
      <value value="500"/>
      <value value="1000"/>
      <value value="2000"/>
      <value value="4000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disagg_rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agg_rate">
      <value value="0.75"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="disagg_rate" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="300"/>
    <metric>soluble_conc</metric>
    <enumeratedValueSet variable="weight0">
      <value value="15"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportion1">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="weight1">
      <value value="64"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="linkwise?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temperature">
      <value value="70"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calculate_dynamically?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stickiness">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disagg_rate">
      <value value="0.005"/>
      <value value="0.01"/>
      <value value="0.015"/>
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population">
      <value value="4000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agg_rate">
      <value value="0.75"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="stickiness" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="250"/>
    <metric>soluble_conc</metric>
    <enumeratedValueSet variable="weights">
      <value value="&quot;[12 24 36 48 60 72 84 100]&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="linkwise?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temperature">
      <value value="80"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportions">
      <value value="&quot;[0.1 0.15 0.3 0.2 0.1 0.05 0.05 0.05]&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calculate_dynamically?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stickiness">
      <value value="0.8"/>
      <value value="1.3"/>
      <value value="1.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="disagg_rate">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="population">
      <value value="200"/>
      <value value="500"/>
      <value value="1000"/>
      <value value="2000"/>
      <value value="4000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="agg_rate">
      <value value="0.75"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
