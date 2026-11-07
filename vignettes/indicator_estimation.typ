// Simple numbering for non-book documents
#let equation-numbering = "(1)"
#let callout-numbering = "1"
#let subfloat-numbering(n-super, subfloat-idx) = {
  numbering("1a", n-super, subfloat-idx)
}

// Theorem configuration for theorion
// Simple numbering for non-book documents (no heading inheritance)
#let theorem-inherited-levels = 0

// Theorem numbering format (can be overridden by extensions for appendix support)
// This function returns the numbering pattern to use
#let theorem-numbering(loc) = "1.1"

// Default theorem render function
#let theorem-render(prefix: none, title: "", full-title: auto, body) = {
  if full-title != "" and full-title != auto and full-title != none {
    strong[#full-title.]
    h(0.5em)
  }
  body
}
// Some definitions presupposed by pandoc's typst output.
#let content-to-string(content) = {
  if content.has("text") {
    content.text
  } else if content.has("children") {
    content.children.map(content-to-string).join("")
  } else if content.has("body") {
    content-to-string(content.body)
  } else if content == [ ] {
    " "
  }
}

#let horizontalrule = line(start: (25%,0%), end: (75%,0%))

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms.item: it => block(breakable: false)[
  #text(weight: "bold")[#it.term]
  #block(inset: (left: 1.5em, top: -0.4em))[#it.description]
]

// Some quarto-specific definitions.

#show raw.where(block: true): set block(
    fill: luma(230),
    width: 100%,
    inset: 8pt,
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let fields = old_block.fields()
  let _ = fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.abs
  }
  block.with(..fields)(new_content)
}

#let empty(v) = {
  if type(v) == str {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == content {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

// Subfloats
// This is a technique that we adapted from https://github.com/tingerrr/subpar/
#let quartosubfloatcounter = counter("quartosubfloatcounter")

#let quarto_super(
  kind: str,
  caption: none,
  label: none,
  supplement: str,
  position: none,
  subcapnumbering: "(a)",
  body,
) = {
  context {
    let figcounter = counter(figure.where(kind: kind))
    let n-super = figcounter.get().first() + 1
    set figure.caption(position: position)
    [#figure(
      kind: kind,
      supplement: supplement,
      caption: caption,
      {
        show figure.where(kind: kind): set figure(numbering: _ => {
          let subfloat-idx = quartosubfloatcounter.get().first() + 1
          subfloat-numbering(n-super, subfloat-idx)
        })
        show figure.where(kind: kind): set figure.caption(position: position)

        show figure: it => {
          let num = numbering(subcapnumbering, n-super, quartosubfloatcounter.get().first() + 1)
          show figure.caption: it => block({
            num.slice(2) // I don't understand why the numbering contains output that it really shouldn't, but this fixes it shrug?
            [ ]
            it.body
          })

          quartosubfloatcounter.step()
          it
          counter(figure.where(kind: it.kind)).update(n => n - 1)
        }

        quartosubfloatcounter.update(0)
        body
      }
    )#label]
  }
}

// callout rendering
// this is a figure show rule because callouts are crossreferenceable
#show figure: it => {
  if type(it.kind) != str {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let children = old_title_block.body.body.children
  let old_title = if children.len() == 1 {
    children.at(0)  // no icon: title at index 0
  } else {
    children.at(1)  // with icon: title at index 1
  }

  // TODO use custom separator if available
  // Use the figure's counter display which handles chapter-based numbering
  // (when numbering is a function that includes the heading counter)
  let callout_num = it.counter.display(it.numbering)
  let new_title = if empty(old_title) {
    [#kind #callout_num]
  } else {
    [#kind #callout_num: #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block,
    block_with_new_content(
      old_title_block.body,
      if children.len() == 1 {
        new_title  // no icon: just the title
      } else {
        children.at(0) + new_title  // with icon: preserve icon block + new title
      }))

  align(left, block_with_new_content(old_callout,
    block(below: 0pt, new_title_block) +
    old_callout.body.children.at(1)))
}

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black, body_background_color: white) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color,
        width: 100%,
        inset: 8pt)[#if icon != none [#text(icon_color, weight: 900)[#icon] ]#title]) +
      if(body != []){
        block(
          inset: 1pt, 
          width: 100%, 
          block(fill: body_background_color, width: 100%, inset: 8pt, body))
      }
    )
}


// syntax highlighting functions from skylighting:
/* Function definitions for syntax highlighting generated by skylighting: */
#let EndLine() = raw("\n")
#let Skylighting(fill: none, number: false, start: 1, sourcelines) = {
   let blocks = []
   let lnum = start - 1
   let bgcolor = rgb("#f1f3f5")
   for ln in sourcelines {
     if number {
       lnum = lnum + 1
       blocks = blocks + box(width: if start + sourcelines.len() > 999 { 30pt } else { 24pt }, text(fill: rgb("#aaaaaa"), [ #lnum ]))
     }
     blocks = blocks + ln + EndLine()
   }
   block(fill: bgcolor, width: 100%, inset: 8pt, radius: 2pt, blocks)
}
#let AlertTok(s) = text(fill: rgb("#ad0000"),raw(s))
#let AnnotationTok(s) = text(fill: rgb("#5e5e5e"),raw(s))
#let AttributeTok(s) = text(fill: rgb("#657422"),raw(s))
#let BaseNTok(s) = text(fill: rgb("#ad0000"),raw(s))
#let BuiltInTok(s) = text(fill: rgb("#003b4f"),raw(s))
#let CharTok(s) = text(fill: rgb("#20794d"),raw(s))
#let CommentTok(s) = text(fill: rgb("#5e5e5e"),raw(s))
#let CommentVarTok(s) = text(style: "italic",fill: rgb("#5e5e5e"),raw(s))
#let ConstantTok(s) = text(fill: rgb("#8f5902"),raw(s))
#let ControlFlowTok(s) = text(weight: "bold",fill: rgb("#003b4f"),raw(s))
#let DataTypeTok(s) = text(fill: rgb("#ad0000"),raw(s))
#let DecValTok(s) = text(fill: rgb("#ad0000"),raw(s))
#let DocumentationTok(s) = text(style: "italic",fill: rgb("#5e5e5e"),raw(s))
#let ErrorTok(s) = text(fill: rgb("#ad0000"),raw(s))
#let ExtensionTok(s) = text(fill: rgb("#003b4f"),raw(s))
#let FloatTok(s) = text(fill: rgb("#ad0000"),raw(s))
#let FunctionTok(s) = text(fill: rgb("#4758ab"),raw(s))
#let ImportTok(s) = text(fill: rgb("#00769e"),raw(s))
#let InformationTok(s) = text(fill: rgb("#5e5e5e"),raw(s))
#let KeywordTok(s) = text(weight: "bold",fill: rgb("#003b4f"),raw(s))
#let NormalTok(s) = text(fill: rgb("#003b4f"),raw(s))
#let OperatorTok(s) = text(fill: rgb("#5e5e5e"),raw(s))
#let OtherTok(s) = text(fill: rgb("#003b4f"),raw(s))
#let PreprocessorTok(s) = text(fill: rgb("#ad0000"),raw(s))
#let RegionMarkerTok(s) = text(fill: rgb("#003b4f"),raw(s))
#let SpecialCharTok(s) = text(fill: rgb("#5e5e5e"),raw(s))
#let SpecialStringTok(s) = text(fill: rgb("#20794d"),raw(s))
#let StringTok(s) = text(fill: rgb("#20794d"),raw(s))
#let VariableTok(s) = text(fill: rgb("#111111"),raw(s))
#let VerbatimStringTok(s) = text(fill: rgb("#20794d"),raw(s))
#let WarningTok(s) = text(style: "italic",fill: rgb("#5e5e5e"),raw(s))



#let article(
  title: none,
  subtitle: none,
  authors: none,
  keywords: (),
  date: none,
  abstract-title: none,
  abstract: none,
  thanks: none,
  cols: 1,
  lang: "en",
  region: "US",
  font: none,
  fontsize: 11pt,
  title-size: 1.5em,
  subtitle-size: 1.25em,
  heading-family: none,
  heading-weight: "bold",
  heading-style: "normal",
  heading-color: black,
  heading-line-height: 0.65em,
  mathfont: none,
  codefont: none,
  linestretch: 1,
  sectionnumbering: none,
  linkcolor: none,
  citecolor: none,
  filecolor: none,
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  // Set document metadata for PDF accessibility
  set document(title: title, keywords: keywords)
  set document(
    author: authors.map(author => content-to-string(author.name)).join(", ", last: " & "),
  ) if authors != none and authors != ()
  set par(
    justify: true,
    leading: linestretch * 0.65em
  )
  set text(lang: lang,
           region: region,
           size: fontsize)
  set text(font: font) if font != none
  show math.equation: set text(font: mathfont) if mathfont != none
  show raw: set text(font: codefont) if codefont != none

  set heading(numbering: sectionnumbering)

  show link: set text(fill: rgb(content-to-string(linkcolor))) if linkcolor != none
  show ref: set text(fill: rgb(content-to-string(citecolor))) if citecolor != none
  show link: this => {
    if filecolor != none and type(this.dest) == label {
      text(this, fill: rgb(content-to-string(filecolor)))
    } else {
      text(this)
    }
   }

  let has-title-block = title != none or (authors != none and authors != ()) or date != none or abstract != none
  if has-title-block {
    place(
      top,
      float: true,
      scope: "parent",
      clearance: 4mm,
      block(below: 1em, width: 100%)[

        #if title != none {
          align(center, block(inset: 2em)[
            #set par(leading: heading-line-height) if heading-line-height != none
            #set text(font: heading-family) if heading-family != none
            #set text(weight: heading-weight)
            #set text(style: heading-style) if heading-style != "normal"
            #set text(fill: heading-color) if heading-color != black

            #text(size: title-size)[#title #if thanks != none {
              footnote(thanks, numbering: "*")
              counter(footnote).update(n => n - 1)
            }]
            #(if subtitle != none {
              parbreak()
              text(size: subtitle-size)[#subtitle]
            })
          ])
        }

        #if authors != none and authors != () {
          let count = authors.len()
          let ncols = calc.min(count, 3)
          grid(
            columns: (1fr,) * ncols,
            row-gutter: 1.5em,
            ..authors.map(author =>
                align(center)[
                  #author.name \
                  #author.affiliation \
                  #author.email
                ]
            )
          )
        }

        #if date != none {
          align(center)[#block(inset: 1em)[
            #date
          ]]
        }

        #if abstract != none {
          block(inset: 2em)[
          #text(weight: "semibold")[#abstract-title] #h(1em) #abstract
          ]
        }
      ]
    )
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
    #outline(
      title: toc_title,
      depth: toc_depth,
      indent: toc_indent
    );
    ]
  }

  doc
}

#set table(
  inset: 6pt,
  stroke: none
)
#import "@preview/fontawesome:0.5.0": *
#let brand-color = (:)
#let brand-color-background = (:)
#let brand-logo = (:)

#set page(
  paper: "us-letter",
  margin: (x: 1.25in, y: 1.25in),
  numbering: "1",
  columns: 1,
)

#show: doc => article(
  title: [Calculating Scaled Ecosystem Condition Indicators from Species Occupancy],
  authors: (
    ( name: [Kwaku Peprah Adjei],
      affiliation: [],
      email: [] ),
    ),
  date: [2026-11-07],
  toc_title: [Table of contents],
  toc_depth: 3,
  doc,
)

= Introduction and Overview
<introduction-and-overview>
== Introduction
<introduction>
Ecosystem condition describes the quality of an ecosystem in terms of its physical, chemical, structural, compositional, and functional characteristics, and its capacity to maintain ecological processes and support biodiversity over time @maes2018@ipbes2019. Assessment of ecosystem condition generally requires the use of indicators, which are measurable variables that summarize attributes of ecosystem structure, composition, or functioning @noss1990. Such indicators provide a practical means of monitoring ecosystems and quantifying changes resulting from natural or anthropogenic pressures.

Breeding bird communities integrate information across multiple taxa and are widely used as indicators of ecosystem composition @hooper2005. In this workflow, species occupancy predictions (see #link("https://github.com/NINAnor/breedingBirdIndicator_OccupancyModelling/blob/main/vignettes/occupancy_modelling.qmd")[vignette]) are scaled against reference values and combined across species into an indicator of ecosystem condition.

== Reference conditions and scaling
<reference-conditions-and-scaling>
Assessing ecosystem condition requires comparison with a reference condition. Reference conditions represent ecosystem states against which current conditions can be evaluated and may correspond to natural, historical, least-disturbed, or management-defined states @stoddard2006. A reference value is the numerical value of an indicator associated with the reference condition. Current ecosystem condition can then be expressed relative to the reference value, providing a standardized measure of ecological condition.

The indicator described here uses #strong[expert-derived reference conditions for each species]. For every species, experts assigned a reference category that expresses how the species' occupancy under reference conditions compares with its predicted occupancy across the region. Each category is linked to a fixed quantile of the species' predicted occupancy distribution within a nature region. The value of that quantile is the species' reference value (called #NormalTok("x100"); below) in that region. Predicted occupancy at every location is then scaled against its reference value, so that a location with occupancy equal to or above the reference value receives the maximum scaled value of 1.

An earlier approach, which used protected areas as proxies for reference conditions and scaled expected species richness (the sum of occupancy probabilities across species), is still available and is described briefly at the end of this vignette (#ref(<sec-protected>, supplement: [Section])).

== Overview
<overview>
The workflow consists of the following steps:

+ Obtain the expert reference quantile for each species.
+ Load the nature regions used for regional reference values and summaries.
+ Scale the occupancy prediction of every species and location against the species' regional reference value (#NormalTok("x100");).
+ Average the scaled values across species at each location.
+ Summarise the indicator for each region and nationally with bootstrapped beta models.
+ Assess the sensitivity of the results to the reference quantile.
+ Save and plot the results.

The indicator is estimated separately for each group of best models (multi-species and single-species models) and each ecosystem type.

= Expert reference quantiles
<expert-reference-quantiles>
== Reference categories and quantiles
<reference-categories-and-quantiles>
The experts' reference categories are read from the species reference workbook for the ecosystem (for example #NormalTok("mountain_experts_refs.xlsx");) and matched to the species list (#NormalTok("SPECIES_SELECTION.xlsx");) by Norwegian species name. Species without a reference category are excluded from the indicator. The experts were given the chance to look at the average occupancy probability estimated from the model and were asked if the occupancy probability of the species will be higher (#NormalTok("høyere");), lower (#NormalTok("lavere"); ) or the same (#NormalTok("samme");) as is in the model given a perfect detection.

Each category is converted to a quantile with #NormalTok("reference_quantile_lookup()");. The defaults are:

#table(
  columns: 2,
  align: (auto,auto,),
  table.header([Reference category], [Reference quantile],),
  table.hline(),
  [#NormalTok("samme");], [0.75],
  [#NormalTok("høyere");], [0.95],
)
These values are fixed within an analysis. A quantile can be changed by naming the category, for example #NormalTok("reference_quantile_lookup(samme = 0.6)");. All other categories keep their defaults. Quantiles must be strictly between 0 and 1.

#NormalTok("get_reference_quantiles()"); reads the workbooks and returns one reference quantile per species:

#Skylighting(([#NormalTok("reference_data ");#OtherTok("<-");#NormalTok(" ");#FunctionTok("get_reference_quantiles");#NormalTok("(");],
[#NormalTok("  ");#AttributeTok("expert_folder =");#NormalTok(" expertFolder,");],
[#NormalTok("  ");#AttributeTok("ecosystem_type =");#NormalTok(" ecosystemType,");],
[#NormalTok("  ");#AttributeTok("ref_quantiles =");#NormalTok(" ");#FunctionTok("reference_quantile_lookup");#NormalTok("()");],
[#NormalTok(")");],));
The function warns about categories that are not in the lookup and about expert species that are missing from the species list, because these species would otherwise be dropped silently.

== Nature regions
<nature-regions>
Reference values are calculated separately for each nature region, which accounts for natural variation in occupancy across environmental gradients. The regions are derived from the county boundaries, aligned to the model coordinate reference system, and saved as a GeoParquet file:

#Skylighting(([#NormalTok("bio_regions ");#OtherTok("<-");#NormalTok(" ");#FunctionTok("get_nature_regions");#NormalTok("(");],
[#NormalTok("  ");#AttributeTok("input_file =");#NormalTok(" gdb,");],
[#NormalTok("  ");#AttributeTok("modelCrs =");#NormalTok(" modelCrs,");],
[#NormalTok("  ");#AttributeTok("output_file =");#NormalTok(" ");#FunctionTok("file.path");#NormalTok("(");],
[#NormalTok("    dataFolder,");],
[#NormalTok("    ");#StringTok("\"processedData\"");#NormalTok(",");],
[#NormalTok("    ");#StringTok("\"nature_regions\"");#NormalTok(",");],
[#NormalTok("    ");#StringTok("\"nature_regions.parquet\"");],
[#NormalTok("  ),");],
[#NormalTok("  ");#AttributeTok("overwrite =");#NormalTok(" ");#ConstantTok("FALSE");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("plot =");#NormalTok(" ");#ConstantTok("TRUE");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("plot_file =");#NormalTok(" ");#FunctionTok("file.path");#NormalTok("(");],
[#NormalTok("    dataFolder,");],
[#NormalTok("    ");#StringTok("\"processedData\"");#NormalTok(",");],
[#NormalTok("    ");#StringTok("\"nature_regions\"");#NormalTok(",");],
[#NormalTok("    ");#StringTok("\"figures\"");],
[#NormalTok("  )");],
[#NormalTok(")");],));
= Scaling species occupancy
<scaling-species-occupancy>
== Reference value (x100)
<reference-value-x100>
Predictions are the posterior median or mean occupancy probability of each species at every prediction location. Locations are assigned to nature regions once, because all species share the same prediction locations. A location that falls in more than one region is assigned to the first; locations outside all regions are given no region and excluded from the regional summaries.

For species $s$ with reference quantile $q_s$ and region $r$, the reference value is the $q_s$-quantile of predicted occupancy over all locations in the region:

#math.equation(block: true, numbering: equation-numbering, [ $ x_(100\,s\,r) = Q_(q_s) (psi_(s\,i) : i in r)\, $ ])<eq-x100>

where $psi_(s\,i)$ is the predicted mean occupancy (#NormalTok("occupancy_mean");) of species $s$ at location $i$. A reference value is set to missing, and a warning is issued, if it is not positive.

== Scaling
<scaling>
Occupancy is scaled with #NormalTok("ecTools::ec_normalise()"); using $x_(100\,s\,r)$ as the value that corresponds to the reference condition. With the lower bound $x_0 = 0$ and the default linear function, the scaled value is

#math.equation(block: true, numbering: equation-numbering, [ $ z_(s\,i) = min (1 \, #h(0em) max (0 \, #h(0em) frac(psi_(s\,i) - x_0, x_(100\,s\,r) - x_0)))\, $ ])<eq-scaled>

so that $z_(s\,i) = 1$ whenever occupancy is at or above the reference value, and $z_(s\,i)$ decreases linearly to 0 as occupancy decreases to $x_0$. Truncation to $\[0\,1\]$ can be switched off with #NormalTok("truncation = FALSE");.

The scaling is implemented in #NormalTok("scale_species_occupancy_by_region()");, which returns

- an #NormalTok("sf"); object with the scaled value of every species at every location, together with the nature region and geometry, and
- a table of the reference values ($x_100$) for each region and species.

== Mean across species
<mean-across-species>
The location-level indicator is the mean of the scaled species values:

#math.equation(block: true, numbering: equation-numbering, [ $ I_i = 1 / S_i sum_(s = 1)^(S_i) z_(s\,i)\, $ ])<eq-indicator>

where $S_i$ is the number of species with a value at location $i$. Missing species values are ignored by default (#NormalTok("na.rm = TRUE");). Because each species is scaled against its own reference value, species that are naturally common and species that are naturally rare contribute on the same scale.

#NormalTok("species_specific_scaling()"); reads the species predictions (only the occupancy column and geometry are kept to reduce memory use), calls #NormalTok("scale_species_occupancy_by_region()");, and adds the mean across species as #NormalTok("mean_scaled");.

= Regional and national estimates
<regional-and-national-estimates>
== Regional estimates
<regional-estimates>
For each nature region, the scaled location values $I_i$ are summarised with a bootstrap:

+ Draw $n = min\(n_(upright("Samps"))\,N_r\)$ values with replacement from the $N_r$ valid values of the region.
+ Fit an intercept-only beta regression, $I tilde.op upright("Beta")$, with #NormalTok("betareg::betareg()");.
+ Back-transform the intercept to the response scale with the inverse logit.
+ Repeat for #NormalTok("n.boot"); replicates.

The mean, median, 2.5% and 97.5% quantiles, and standard deviation of the replicate estimates are reported. Beta regression requires values strictly between 0 and 1. Values of exactly 0 or 1, which can occur after truncation, are squeezed into the open interval with the transformation $\(y\(n - 1\)+ 0.5\)\/n$ @smithson2006, applied only when such values are present.

== National estimate
<national-estimate>
The national estimate pools the values of all regions and fits an intercept-only beta mixed model with region as a random effect,

#math.equation(block: true, numbering: equation-numbering, [ $ upright("logit")\(mu_(i r)\)= beta_0 + u_r\,#h(2em) u_r tilde.op N\(0\,sigma_u^2\)\, $ ])<eq-mixed>

with #NormalTok("glmmTMB::glmmTMB(y ~ 1 + (1 | region), family = beta_family())");. The bootstrap is stratified: each replicate draws up to $n_(upright("Samps"))$ values with replacement within each region, so that every region keeps its own weight regardless of its size. The national estimate is the back-transformed fixed intercept $upright("logit")^(- 1)\(beta_0\)$, which describes a typical region (random effect equal to zero). It is therefore not the same as the pooled mean of all locations. Replicates for which the model does not fit are dropped, and a warning reports how many.

#block[
#callout(
body: 
[
The bootstrap treats the prediction locations as independent observations. Neighbouring locations are spatially autocorrelated, and the occupancy predictions themselves are uncertain. The bootstrap intervals therefore describe resampling variability of the locations only and understate the total uncertainty of the indicator.

]
, 
title: 
[
Note
]
, 
background_color: 
rgb("#dae6fb")
, 
icon_color: 
rgb("#0758E5")
, 
icon: 
fa-info()
, 
body_background_color: 
white
)
]
= Running the workflow
<running-the-workflow>
The workflow is run by #NormalTok("03_estimating_indicators_using_expert_references.R");. It is selected in the master script with #NormalTok("useExpertData <- TRUE");. The script

+ gets the reference quantiles and the nature regions,
+ finds the prediction folders of the best models (#NormalTok("prediction"); folder, with #NormalTok("prediction_outputs"); as a fallback for older runs),
+ runs #NormalTok("calculate_scaled_ecosystem_richness_with_expert_data()"); for the multi-species and single-species groups, and
+ saves the results and a combined summary across the model groups.

#Skylighting(([#ControlFlowTok("if");#NormalTok(" (useExpertData) {");],
[#NormalTok("  ");#FunctionTok("source");#NormalTok("(here");#SpecialCharTok("::");#FunctionTok("here");#NormalTok("(");],
[#NormalTok("    ");#StringTok("\"scripts\"");#NormalTok(",");],
[#NormalTok("    ");#StringTok("\"02_pipeline\"");#NormalTok(",");],
[#NormalTok("    ");#StringTok("\"02_modelFit\"");#NormalTok(",");],
[#NormalTok("    ");#StringTok("\"03_estimating_indicators_using_expert_references.R\"");],
[#NormalTok("  ))");],
[#NormalTok("} ");#ControlFlowTok("else");#NormalTok(" {");],
[#NormalTok("  ");#FunctionTok("source");#NormalTok("(here");#SpecialCharTok("::");#FunctionTok("here");#NormalTok("(");],
[#NormalTok("    ");#StringTok("\"scripts\"");#NormalTok(",");],
[#NormalTok("    ");#StringTok("\"02_pipeline\"");#NormalTok(",");],
[#NormalTok("    ");#StringTok("\"02_modelFit\"");#NormalTok(",");],
[#NormalTok("    ");#StringTok("\"03_estimating_indicators.R\"");],
[#NormalTok("  ))");],
[#NormalTok("}");],));
The main function can also be called directly:

#Skylighting(([#NormalTok("scaled_indicator ");#OtherTok("<-");#NormalTok(" ");#FunctionTok("calculate_scaled_ecosystem_richness_with_expert_data");#NormalTok("(");],
[#NormalTok("  ");#AttributeTok("input_dir =");#NormalTok(" path_to_predictions,");],
[#NormalTok("  ");#AttributeTok("reference_data =");#NormalTok(" reference_data,");],
[#NormalTok("  ");#AttributeTok("categorisation_shape =");#NormalTok(" bio_regions,");],
[#NormalTok("  ");#AttributeTok("ecosystem_type =");#NormalTok(" ecosystemType,");],
[#NormalTok("  ");#AttributeTok("value_col =");#NormalTok(" ");#StringTok("\"occupancy_mean\"");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("region_col =");#NormalTok(" ");#StringTok("\"nature_region\"");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("output_dir =");#NormalTok(" path_to_indicators,");],
[#NormalTok("  ");#AttributeTok("save_outputs =");#NormalTok(" ");#ConstantTok("TRUE");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("overwrite =");#NormalTok(" ");#ConstantTok("TRUE");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("rep_values_for_each_region =");#NormalTok(" ");#ConstantTok("TRUE");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("n.boot =");#NormalTok(" ");#DecValTok("500");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("boot.seed =");#NormalTok(" ");#DecValTok("123");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("n.Samps =");#NormalTok(" ");#DecValTok("20000");],
[#NormalTok(")");],));
== Outputs
<outputs>
The function returns a list and, if #NormalTok("save_outputs = TRUE");, saves each element in the output folder, with GeoParquet for spatial data and csv for tables.

#table(
  columns: (37.5%, 25%, 37.5%),
  align: (auto,auto,auto,),
  table.header([Element], [File], [Content],),
  table.hline(),
  [#NormalTok("indicator_vals");], [#NormalTok("indicator_vals.parquet");], [Per location: unscaled mean occupancy across species (#NormalTok("indicator_vals");), #NormalTok("scaled_indicator"); (#ref(<eq-indicator>, supplement: [Equation])), region and geometry],
  [#NormalTok("scaled_species");], [#NormalTok("scaled_species.parquet");], [Per location: scaled value of every species (#ref(<eq-scaled>, supplement: [Equation])), #NormalTok("mean_scaled");, region and geometry],
  [#NormalTok("x100_values");], [#NormalTok("x100_values.csv");], [Reference value $x_100$ (#ref(<eq-x100>, supplement: [Equation])) for every region and species],
  [#NormalTok("species_summary");], [#NormalTok("species_summary.csv");], [Mean, median, 2.5% and 97.5% quantiles of the scaled values of every species, by region and nationally],
  [#NormalTok("representative_indicator_vals");], [#NormalTok("representative_indicator_vals.csv");], [Bootstrap mean, median, 2.5% and 97.5% quantiles, and standard deviation of the indicator for each region (beta regression) and #NormalTok("National"); (beta mixed model)],
  [#NormalTok("reference_data");], [#NormalTok("reference_data.csv");], [Reference quantiles of the species used],
  [all of the above], [#NormalTok("all_results.rds");], [Complete output object],
)
When models are compared (#NormalTok("compareMultipleOutputs <- TRUE");), results are saved in #NormalTok("indicators_outputs/compareMultipleOutputs/<ecosystem>/<model group>/");. Combined tables across groups (#NormalTok("indicator_summary.csv");, #NormalTok("x100_values_all_groups.csv");, #NormalTok("species_summary_all_groups.csv");) are saved in the folder above.

== Maps and summary plots
<maps-and-summary-plots>
#NormalTok("01_process_one_path_for_results.R"); reads the saved results and creates maps of the unscaled mean occupancy and the scaled indicator, a plot of the regional and national estimates with their 95% intervals, and a tile plot of the $x_100$ values by species and region.

= Sensitivity to the reference quantile
<sensitivity-to-the-reference-quantile>
The choice of quantile for a reference category is a subjective decision. #NormalTok("reference_quantile_sensitivity()"); repeats the analysis for several quantiles of one category (default #NormalTok("samme");: 0.5, 0.6, 0.7, 0.8 and 0.9; the other category keeps its default). The species predictions are read once, and the same bootstrap seed is used for every quantile so that differences reflect the quantile and not bootstrap noise.

#Skylighting(([#NormalTok("sensitivity_results ");#OtherTok("<-");#NormalTok(" ");#FunctionTok("reference_quantile_sensitivity");#NormalTok("(");],
[#NormalTok("  ");#AttributeTok("quantile_values =");#NormalTok(" ");#FunctionTok("c");#NormalTok("(");#FloatTok("0.5");#NormalTok(", ");#FloatTok("0.6");#NormalTok(", ");#FloatTok("0.7");#NormalTok(", ");#FloatTok("0.8");#NormalTok(", ");#FloatTok("0.9");#NormalTok("),");],
[#NormalTok("  ");#AttributeTok("category =");#NormalTok(" ");#StringTok("\"samme\"");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("input_dir =");#NormalTok(" sens_pred_dir,");],
[#NormalTok("  ");#AttributeTok("expert_folder =");#NormalTok(" expertFolder,");],
[#NormalTok("  ");#AttributeTok("ecosystem_type =");#NormalTok(" ecosystemType,");],
[#NormalTok("  ");#AttributeTok("categorisation_shape =");#NormalTok(" bio_regions,");],
[#NormalTok("  ");#AttributeTok("region_col =");#NormalTok(" ");#StringTok("\"nature_region\"");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("n.boot =");#NormalTok(" ");#DecValTok("100");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("n.Samps =");#NormalTok(" ");#DecValTok("5000");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("output_dir =");#NormalTok(" sens_dir");],
[#NormalTok(")");],
[],
[#NormalTok("sensitivity_plots ");#OtherTok("<-");#NormalTok(" ");#FunctionTok("plot_reference_quantile_sensitivity");#NormalTok("(");],
[#NormalTok("  sensitivity_results,");],
[#NormalTok("  ");#AttributeTok("category =");#NormalTok(" ");#StringTok("\"samme\"");],
[#NormalTok(")");],));
The tables (#NormalTok("sensitivity_x100_values.csv");, #NormalTok("sensitivity_representative_indicator_vals.csv");, #NormalTok("sensitivity_species_summary.csv");) and two figures are saved in #NormalTok("indicators_outputs/sensitivity/<ecosystem>/samme_quantile/");. The figures show the indicator estimate and its interval against the quantile for each region, and $x_100$ against the quantile for each species of the varied category.

In the mountain ecosystem, raising the #NormalTok("samme"); quantile from 0.5 to 0.9 raised the reference values of the affected species and lowered the national indicator from about 0.64 to about 0.54, which is about 0.024 per 0.1 increase in the quantile. The change was similar in all regions, so the ranking of regions was largely unaffected. The change came mostly from species with low occupancy over most of a region, for which the quantile moves from near zero to a moderate value, whereas species whose reference value was already close to 1 were unaffected. These results depend on the species set and predictions of that analysis.

= Interpretation
<interpretation>
The indicator expresses how closely predicted occupancy approaches the expert-defined reference condition, averaged across species:

- values close to $1$ indicate that, on average, species' occupancy is at or above their reference values;
- values close to $0$ indicate low occupancy relative to the reference values.

Because $x_100$ is a quantile of the predictions in the region, the indicator is relative to the occupancy that the model predicts within that region and to the quantile chosen by the experts' categories. A species whose reference value is very low (a rare species in that region) reaches the maximum scaled value almost everywhere, and a species whose reference value is close to 1 contributes little variation. The results should therefore be read together with the $x_100$ values and the sensitivity analysis.

Because the indicator is derived from occupancy probabilities rather than species inventories, it should be interpreted as a model-based measure of relative ecosystem condition rather than a direct census of biodiversity.

= Alternative: protected areas as reference conditions
<sec-protected>
The earlier approach uses protected areas as proxies for reference conditions and scales #strong[expected richness], the sum of species occupancy probabilities, rather than species-specific occupancy @dorazio2006. It is run by #NormalTok("03_estimating_indicators.R"); when #NormalTok("useExpertData <- FALSE");.

Protected areas are loaded with #NormalTok("get_protected_areas()");, and a regular grid of points inside the protected-area polygons is used to extract species predictions (#NormalTok("extract_species_values_at_reference()");). Expected richness at a location is

#Skylighting(([#NormalTok("expected_richness ");#OtherTok("=");#NormalTok(" ");#FunctionTok("rowSums");#NormalTok("(species_prediction_values, ");#AttributeTok("na.rm =");#NormalTok(" ");#ConstantTok("TRUE");#NormalTok(")");],));
The maximum expected richness observed among the reference points in each biogeographic region is the regional reference value $R_(upright(r e f)\,r)$, and scaled richness is calculated as

$ upright("Scaled richness") = 0.6 times frac(R_i - R_min, R_(upright(r e f)\,r) - R_min)\, $

where $R_i$ is the expected richness at location $i$ and $R_min$ is the minimum expected richness across the prediction domain (assumed to be $0$). Scaled richness values are constrained to the interval $\[0\,1\]$. The reference value corresponds to $0.6$ on the scaled indicator, so values greater than $0.6$ indicate conditions better than the reference condition.

#Skylighting(([#NormalTok("scaled_richness ");#OtherTok("<-");#NormalTok(" ");#FunctionTok("calculate_scaled_ecosystem_richness");#NormalTok("(");],
[#NormalTok("  ");#AttributeTok("input_dir =");#NormalTok(" path_to_predictions,");],
[#NormalTok("  ");#AttributeTok("reference_shape =");#NormalTok(" protected_regions,");],
[#NormalTok("  ");#AttributeTok("categorisation_shape =");#NormalTok(" bio_regions,");],
[#NormalTok("  ");#AttributeTok("ecosystem_type =");#NormalTok(" ecosystemType,");],
[#NormalTok("  ");#AttributeTok("value_col =");#NormalTok(" indicator_col,");],
[#NormalTok("  ");#AttributeTok("point_spacing =");#NormalTok(" aggregation_spacing,");],
[#NormalTok("  ");#AttributeTok("forest_only_strict =");#NormalTok(" ecosystemType ");#SpecialCharTok("==");#NormalTok(" ");#StringTok("\"Forest\"");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("na.rm =");#NormalTok(" ");#ConstantTok("TRUE");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("output_dir =");#NormalTok(" path_to_indicators,");],
[#NormalTok("  ");#AttributeTok("save_outputs =");#NormalTok(" ");#ConstantTok("TRUE");#NormalTok(",");],
[#NormalTok("  ");#AttributeTok("overwrite =");#NormalTok(" ");#ConstantTok("FALSE");],
[#NormalTok(")");],));
For forest ecosystems, the reference condition can be restricted to protected forest ecosystems. When #NormalTok("forest_only_strict = TRUE");, only protected areas classified as national parks or nature reserves and designated as forest protection areas are retained.

#bibliography(("references.bib"))

