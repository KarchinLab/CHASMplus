---
#
# Use the widgets beneath and the content will be
# inserted automagically in the webpage. To make
# this work, you have to use › layout: frontpage
#
layout: frontpage
#header:
#  image_fullwidth: chasmplus.png
#
# Use the call for action to show a button on the frontpage
#
# To make internal links, just use a permalink like this
# url: /getting-started/
#
# To style the button in different colors, use no value
# to use the main color or success, alert or secondary.
# To change colors see sass/_01_settings_colors.scss
#
permalink: /index.html
---
<div style="text-align:center"><h1>Explore driver mutations</h1></div>

<div id="tabulator-controls" class="table-controls">
  <div class="medium-8" style="display:block;margin-left:auto;margin-right:auto;">
    <input name="name" id="geneSearch" type="text" placeholder="Search by gene name" class="form-control">
  </div>
</div>

<div id="bar-chart" class="medium-12 columns"> </div>
<div id="pie-chart" class="medium-6 columns"> </div>
<div id="violin-chart" class="medium-6 columns"> </div>

<div id="tabulator-controls" class="table-controls">
  <div class="medium-3 columns">
    <div style="float:left;">Analysis:</div>
    <select name="ctype" id="analysisDropdown" placeholder="Select an analysis" selected="pan-cancer">
      <option value="pan-cancer">pan-cancer</option>
      <option value="cancer type-specific">cancer type-specific</option>
    </select>
  </div>
  <div class="medium-6 columns" style="display:block;margin-left:auto;margin-right:auto;">
    Columns:<br>
    <div style="white-space: nowrap;display:inline"><input type="checkbox" name="Hugo_Symbol" checked> Gene</div>
    <div style="white-space: nowrap;display:inline"><input type="checkbox" name="Transcript_ID"> Transcript</div>
    <div style="white-space: nowrap;display:inline"><input type="checkbox" name="HGVSp_Short" checked> Mutation</div>
    <div style="white-space: nowrap;display:inline"><input type="checkbox" name="gwCHASMplus score" checked> Score</div>
    <div style="white-space: nowrap;display:inline"><input type="checkbox" name="cancer type" checked> Cancer type</div>
    <div style="white-space: nowrap;display:inline"><input type="checkbox" name="number of mutations" checked> Number of mutations</div>
    <div style="white-space: nowrap;display:inline"><input type="checkbox" name="frequency category" checked> Frequency</div>
    <div style="white-space: nowrap;display:inline"><input type="checkbox" name="url" checked> Detailed information</div>
  </div>
  <div class="medium-3 columns">
    Download:<br>
    <button type="button" name="csv-download">CSV</button>
    <button type="button" name="xlsx-download">Excel</button>
  </div>
</div>
<div id="example-table" class="medium-12 columns"></div>
