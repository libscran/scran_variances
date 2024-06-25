<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.8">
  <compound kind="file">
    <name>choose_highly_variable_genes.hpp</name>
    <path>scran/</path>
    <filename>choose__highly__variable__genes_8hpp.html</filename>
    <namespace>scran</namespace>
    <namespace>scran::choose_highly_variable_genes</namespace>
  </compound>
  <compound kind="file">
    <name>fit_variance_trend.hpp</name>
    <path>scran/</path>
    <filename>fit__variance__trend_8hpp.html</filename>
    <class kind="struct">scran::fit_variance_trend::Options</class>
    <class kind="struct">scran::fit_variance_trend::Workspace</class>
    <class kind="struct">scran::fit_variance_trend::Results</class>
    <namespace>scran</namespace>
    <namespace>scran::fit_variance_trend</namespace>
  </compound>
  <compound kind="file">
    <name>model_gene_variances.hpp</name>
    <path>scran/</path>
    <filename>model__gene__variances_8hpp.html</filename>
    <includes id="fit__variance__trend_8hpp" name="fit_variance_trend.hpp" local="yes" import="no" module="no" objc="no">fit_variance_trend.hpp</includes>
    <class kind="struct">scran::model_gene_variances::Options</class>
    <class kind="struct">scran::model_gene_variances::Results</class>
    <class kind="struct">scran::model_gene_variances::BlockResults</class>
    <namespace>scran</namespace>
    <namespace>scran::model_gene_variances</namespace>
  </compound>
  <compound kind="file">
    <name>scran.hpp</name>
    <path>scran/</path>
    <filename>scran_8hpp.html</filename>
    <namespace>scran</namespace>
  </compound>
  <compound kind="struct">
    <name>scran::model_gene_variances::BlockResults</name>
    <filename>structscran_1_1model__gene__variances_1_1BlockResults.html</filename>
    <templarg>typename Stat_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Results&lt; Stat_ &gt; &gt;</type>
      <name>per_block</name>
      <anchorfile>structscran_1_1model__gene__variances_1_1BlockResults.html</anchorfile>
      <anchor>a8d112d3b54a5dc43f6daa1eb64d89905</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Results&lt; Stat_ &gt;</type>
      <name>average</name>
      <anchorfile>structscran_1_1model__gene__variances_1_1BlockResults.html</anchorfile>
      <anchor>a7849327e43c2ba92301da3e43b6a9db7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran::fit_variance_trend::Options</name>
    <filename>structscran_1_1fit__variance__trend_1_1Options.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>minimum_mean</name>
      <anchorfile>structscran_1_1fit__variance__trend_1_1Options.html</anchorfile>
      <anchor>a594d757567e8546ef4fe7f28fda9b726</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>mean_filter</name>
      <anchorfile>structscran_1_1fit__variance__trend_1_1Options.html</anchorfile>
      <anchor>a7c75573f56a2506d3b0f47cbb1465d22</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>transform</name>
      <anchorfile>structscran_1_1fit__variance__trend_1_1Options.html</anchorfile>
      <anchor>a77a862e1aba3394c8a6c896dfef5cfc3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>span</name>
      <anchorfile>structscran_1_1fit__variance__trend_1_1Options.html</anchorfile>
      <anchor>ade37259a86a113b6dbed01cece15ac64</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>use_minimum_width</name>
      <anchorfile>structscran_1_1fit__variance__trend_1_1Options.html</anchorfile>
      <anchor>ab958e99903b5e586a869280e8fc0bc17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>minimum_width</name>
      <anchorfile>structscran_1_1fit__variance__trend_1_1Options.html</anchorfile>
      <anchor>a01ee328deb637aab3656e3d77dcdd150</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>minimum_window_count</name>
      <anchorfile>structscran_1_1fit__variance__trend_1_1Options.html</anchorfile>
      <anchor>a4e0a6b25f6beade1362c7d5257318947</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran_1_1fit__variance__trend_1_1Options.html</anchorfile>
      <anchor>a0d05231069087607ea1b8c68d4e934e6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran::model_gene_variances::Options</name>
    <filename>structscran_1_1model__gene__variances_1_1Options.html</filename>
    <member kind="variable">
      <type>fit_variance_trend::Options</type>
      <name>fit_variance_trend_options</name>
      <anchorfile>structscran_1_1model__gene__variances_1_1Options.html</anchorfile>
      <anchor>afbf016616201702fc924fd4473fd0d03</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>block_weights::Policy</type>
      <name>block_weight_policy</name>
      <anchorfile>structscran_1_1model__gene__variances_1_1Options.html</anchorfile>
      <anchor>a2e846662cc00387e3c16751dd20589fd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>block_weights::VariableParameters</type>
      <name>variable_block_weight_parameters</name>
      <anchorfile>structscran_1_1model__gene__variances_1_1Options.html</anchorfile>
      <anchor>ad181c41f6e57d0a2f24a5aaf4b4d88d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_average</name>
      <anchorfile>structscran_1_1model__gene__variances_1_1Options.html</anchorfile>
      <anchor>a3fbb5a99b57e24934b6bfcb38014af20</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran_1_1model__gene__variances_1_1Options.html</anchorfile>
      <anchor>a41e193f118e86a571f87528fcfe395e9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran::fit_variance_trend::Results</name>
    <filename>structscran_1_1fit__variance__trend_1_1Results.html</filename>
    <templarg>typename Float_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Float_ &gt;</type>
      <name>fitted</name>
      <anchorfile>structscran_1_1fit__variance__trend_1_1Results.html</anchorfile>
      <anchor>ae3d228637819c6dfdb83d4c63672782b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Float_ &gt;</type>
      <name>residuals</name>
      <anchorfile>structscran_1_1fit__variance__trend_1_1Results.html</anchorfile>
      <anchor>a5f8ad5835bae3bbdd334862ace478ca6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran::model_gene_variances::Results</name>
    <filename>structscran_1_1model__gene__variances_1_1Results.html</filename>
    <templarg>typename Stat_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>means</name>
      <anchorfile>structscran_1_1model__gene__variances_1_1Results.html</anchorfile>
      <anchor>afd573efb30de0e1435a92b7eb1957bcb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>variances</name>
      <anchorfile>structscran_1_1model__gene__variances_1_1Results.html</anchorfile>
      <anchor>a780c73b9c1432a94eba96422405909f1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>fitted</name>
      <anchorfile>structscran_1_1model__gene__variances_1_1Results.html</anchorfile>
      <anchor>a8a64e1d52cec567879579190be5c2b7e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>residuals</name>
      <anchorfile>structscran_1_1model__gene__variances_1_1Results.html</anchorfile>
      <anchor>acda4e88119792728fa27341982798e09</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran::fit_variance_trend::Workspace</name>
    <filename>structscran_1_1fit__variance__trend_1_1Workspace.html</filename>
    <templarg>typename Float_</templarg>
  </compound>
  <compound kind="namespace">
    <name>scran</name>
    <filename>namespacescran.html</filename>
    <namespace>scran::choose_highly_variable_genes</namespace>
    <namespace>scran::fit_variance_trend</namespace>
    <namespace>scran::model_gene_variances</namespace>
  </compound>
  <compound kind="namespace">
    <name>scran::choose_highly_variable_genes</name>
    <filename>namespacescran_1_1choose__highly__variable__genes.html</filename>
    <member kind="function">
      <type>void</type>
      <name>compute</name>
      <anchorfile>namespacescran_1_1choose__highly__variable__genes.html</anchorfile>
      <anchor>a220eab651cc929095ea3bbb52aea03da</anchor>
      <arglist>(size_t n, const Stat_ *statistic, Bool_ *output, const Options &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Bool_ &gt;</type>
      <name>compute</name>
      <anchorfile>namespacescran_1_1choose__highly__variable__genes.html</anchorfile>
      <anchor>a2e18e3a5d19d6c91d23cbc2a4c5633ee</anchor>
      <arglist>(size_t n, const Stat_ *statistic, const Options &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Index_ &gt;</type>
      <name>compute_index</name>
      <anchorfile>namespacescran_1_1choose__highly__variable__genes.html</anchorfile>
      <anchor>a9a0d47de0eb9f7d5623a217611c9615b</anchor>
      <arglist>(Index_ n, const Stat_ *statistic, const Options &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>scran::fit_variance_trend</name>
    <filename>namespacescran_1_1fit__variance__trend.html</filename>
    <class kind="struct">scran::fit_variance_trend::Options</class>
    <class kind="struct">scran::fit_variance_trend::Results</class>
    <class kind="struct">scran::fit_variance_trend::Workspace</class>
    <member kind="function">
      <type>void</type>
      <name>compute</name>
      <anchorfile>namespacescran_1_1fit__variance__trend.html</anchorfile>
      <anchor>ab8296a4054e4eb94a44864c1ad47d53a</anchor>
      <arglist>(size_t n, const Float_ *mean, const Float_ *variance, Float_ *fitted, Float_ *residuals, Workspace&lt; Float_ &gt; &amp;workspace, const Options &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>Results&lt; Float_ &gt;</type>
      <name>compute</name>
      <anchorfile>namespacescran_1_1fit__variance__trend.html</anchorfile>
      <anchor>aa4769d5f59810f4a8b7b8c13a84675ca</anchor>
      <arglist>(size_t n, const Float_ *mean, const Float_ *variance, const Options &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>scran::model_gene_variances</name>
    <filename>namespacescran_1_1model__gene__variances.html</filename>
    <class kind="struct">scran::model_gene_variances::BlockResults</class>
    <class kind="struct">scran::model_gene_variances::Options</class>
    <class kind="struct">scran::model_gene_variances::Results</class>
    <member kind="function">
      <type>void</type>
      <name>compute_blocked</name>
      <anchorfile>namespacescran_1_1model__gene__variances.html</anchorfile>
      <anchor>ad6afe0443f263323ace77cc0d44e1877</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; *mat, const Block_ *block, std::vector&lt; Stat_ * &gt; means, std::vector&lt; Stat_ * &gt; variances, std::vector&lt; Stat_ * &gt; fitted, std::vector&lt; Stat_ * &gt; residuals, Stat_ *ave_means, Stat_ *ave_variances, Stat_ *ave_fitted, Stat_ *ave_residuals, const Options &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_blocked</name>
      <anchorfile>namespacescran_1_1model__gene__variances.html</anchorfile>
      <anchor>a2d39f7cfaaae62be97f2c548498a8b57</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; *mat, const Block_ *block, std::vector&lt; Stat_ * &gt; means, std::vector&lt; Stat_ * &gt; variances, std::vector&lt; Stat_ * &gt; fitted, std::vector&lt; Stat_ * &gt; residuals, const Options &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute</name>
      <anchorfile>namespacescran_1_1model__gene__variances.html</anchorfile>
      <anchor>a341ce22fc3d3aae1752e13b9bb2b304b</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; *mat, Stat_ *means, Stat_ *variances, Stat_ *fitted, Stat_ *residuals, const Options &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>Results&lt; Stat_ &gt;</type>
      <name>compute</name>
      <anchorfile>namespacescran_1_1model__gene__variances.html</anchorfile>
      <anchor>aba82ca062767128111bd1b2b530fafed</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; *mat, const Options &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>BlockResults&lt; Stat_ &gt;</type>
      <name>compute_blocked</name>
      <anchorfile>namespacescran_1_1model__gene__variances.html</anchorfile>
      <anchor>a32290dcde184a4a253b020809532fd35</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; *mat, const Block_ *block, const Options &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Model per-gene variance in expression</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
