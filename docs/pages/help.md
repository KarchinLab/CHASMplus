---
layout: page
title: "Help"
subheadline: ""
teaser: ""
permalink: /help/
---

This mutation portal provides an interactive interface to predictions by CHASMplus.
CHASMplus is a machine learning algorithm that discriminates somatic missense mutations as
cancer drivers. We summarize results with interactive widgets, provide search capability and an interactive table of results. Visualizations include:

* Mutations on protein structure
* Lollipop diagrams
* Number of driver mutations per gene
* Breakdown of mutation frequency
* Distribution of CHASMplus scores

## Search by gene

Enter in the HUGO gene symbol at the top of the page to search for all predicted
driver mutations for a particular gene.

![Figure 1]({{ site.url }}{{ site.baseurl }}/images/gene_search.png)

## Summary plots

All searches will dynamically change the summary plots.

## Table

The table provides a listing of all predicted driver mutations meeting your
search criteria.

![Figure 2]({{ site.url }}{{ site.baseurl }}/images/table.png)

**Additional variant information:**

By clicking on the "More details" link for each variant (see Marker \#8), CRAVAT
provides more information on the variant, including a lollipop plot, protein structure
viewer and gene network graph. You can read help about using these widgets from
the [CRAVAT](http://www.cravat.us/CRAVAT/help.jsp) and [MuPIT](http://mupit.icm.jhu.edu/MuPIT_Interactive/Help.html) help pages.

![Figure 3]({{ site.url }}{{ site.baseurl }}/images/single_variant_page.png)

**Toggle columns:**

You can toggle the display of columns by checking/unchecking boxes (see Marker \#7).

**Search capabilities:**

You may enter multiple combinations of searches by using the advanced searches
available on the table.

1. Filter by cancer type
2. Filter by gene name (exact match)
3. Filter by mutation (fuzzy match)
4. Filter by mutations with at least a certain CHASMplus score
5. Filter by cancer type (in table header)
6. Filter by the mutation frequency of the driver mutation (common, intermediate, or rare)

**Download data:**

You can download data in either excel format (.xlsx file) or as a csv file (see Marker \#9).
