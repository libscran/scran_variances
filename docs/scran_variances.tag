<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.12.0">
  <compound kind="file">
    <name>choose_highly_variable_genes.hpp</name>
    <path>scran_variances/</path>
    <filename>choose__highly__variable__genes_8hpp.html</filename>
    <class kind="struct">scran_variances::ChooseHighlyVariableGenesOptions</class>
    <namespace>scran_variances</namespace>
  </compound>
  <compound kind="file">
    <name>fit_variance_trend.hpp</name>
    <path>scran_variances/</path>
    <filename>fit__variance__trend_8hpp.html</filename>
    <class kind="struct">scran_variances::FitVarianceTrendOptions</class>
    <class kind="struct">scran_variances::FitVarianceTrendWorkspace</class>
    <class kind="struct">scran_variances::FitVarianceTrendResults</class>
    <namespace>scran_variances</namespace>
  </compound>
  <compound kind="file">
    <name>model_gene_variances.hpp</name>
    <path>scran_variances/</path>
    <filename>model__gene__variances_8hpp.html</filename>
    <includes id="fit__variance__trend_8hpp" name="fit_variance_trend.hpp" local="yes" import="no" module="no" objc="no">fit_variance_trend.hpp</includes>
    <class kind="struct">scran_variances::ModelGeneVariancesOptions</class>
    <class kind="struct">scran_variances::ModelGeneVariancesBuffers</class>
    <class kind="struct">scran_variances::ModelGeneVariancesResults</class>
    <class kind="struct">scran_variances::ModelGeneVariancesBlockedBuffers</class>
    <class kind="struct">scran_variances::ModelGeneVariancesBlockedResults</class>
    <namespace>scran_variances</namespace>
  </compound>
  <compound kind="file">
    <name>scran_variances.hpp</name>
    <path>scran_variances/</path>
    <filename>scran__variances_8hpp.html</filename>
    <includes id="fit__variance__trend_8hpp" name="fit_variance_trend.hpp" local="yes" import="no" module="no" objc="no">fit_variance_trend.hpp</includes>
    <includes id="model__gene__variances_8hpp" name="model_gene_variances.hpp" local="yes" import="no" module="no" objc="no">model_gene_variances.hpp</includes>
    <includes id="choose__highly__variable__genes_8hpp" name="choose_highly_variable_genes.hpp" local="yes" import="no" module="no" objc="no">choose_highly_variable_genes.hpp</includes>
    <namespace>scran_variances</namespace>
  </compound>
  <compound kind="struct">
    <name>scran_variances::ChooseHighlyVariableGenesOptions</name>
    <filename>structscran__variances_1_1ChooseHighlyVariableGenesOptions.html</filename>
    <member kind="variable">
      <type>std::size_t</type>
      <name>top</name>
      <anchorfile>structscran__variances_1_1ChooseHighlyVariableGenesOptions.html</anchorfile>
      <anchor>a207b3e1474c88a93e9c9511e3757c43a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>larger</name>
      <anchorfile>structscran__variances_1_1ChooseHighlyVariableGenesOptions.html</anchorfile>
      <anchor>a7ac6ad16838227b7ea5f2266b09c8d37</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>use_bound</name>
      <anchorfile>structscran__variances_1_1ChooseHighlyVariableGenesOptions.html</anchorfile>
      <anchor>a2990bc14c7a1912b0a2bf9d0670be7fb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>bound</name>
      <anchorfile>structscran__variances_1_1ChooseHighlyVariableGenesOptions.html</anchorfile>
      <anchor>a2d75802080acf97aae3f96e4e3e1f2b7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>keep_ties</name>
      <anchorfile>structscran__variances_1_1ChooseHighlyVariableGenesOptions.html</anchorfile>
      <anchor>a452e64f05baa4e9fe2fc35e97510bece</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_variances::FitVarianceTrendOptions</name>
    <filename>structscran__variances_1_1FitVarianceTrendOptions.html</filename>
    <member kind="variable">
      <type>double</type>
      <name>minimum_mean</name>
      <anchorfile>structscran__variances_1_1FitVarianceTrendOptions.html</anchorfile>
      <anchor>aeb49baa7624a6aa1df7099761aa5fe53</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>mean_filter</name>
      <anchorfile>structscran__variances_1_1FitVarianceTrendOptions.html</anchorfile>
      <anchor>acaaa268e4ebcdfcc105fb05cc99fa91c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>transform</name>
      <anchorfile>structscran__variances_1_1FitVarianceTrendOptions.html</anchorfile>
      <anchor>aa1e9b5b27d936b35e9cf18788c66792b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>span</name>
      <anchorfile>structscran__variances_1_1FitVarianceTrendOptions.html</anchorfile>
      <anchor>a5b9937546b3f82e869d0a0563a624fab</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>use_minimum_width</name>
      <anchorfile>structscran__variances_1_1FitVarianceTrendOptions.html</anchorfile>
      <anchor>a563b78705db208dce5727370357e0070</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>minimum_width</name>
      <anchorfile>structscran__variances_1_1FitVarianceTrendOptions.html</anchorfile>
      <anchor>a372f10d4bc53cd8e94297bf50328c3b3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>minimum_window_count</name>
      <anchorfile>structscran__variances_1_1FitVarianceTrendOptions.html</anchorfile>
      <anchor>aaed78ddd1cc4ce209f7b318c67bf22c2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran__variances_1_1FitVarianceTrendOptions.html</anchorfile>
      <anchor>a0811d658119b078fe61a0d65ea367586</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_variances::FitVarianceTrendResults</name>
    <filename>structscran__variances_1_1FitVarianceTrendResults.html</filename>
    <templarg>typename Float_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Float_ &gt;</type>
      <name>fitted</name>
      <anchorfile>structscran__variances_1_1FitVarianceTrendResults.html</anchorfile>
      <anchor>ab3ffa9045cb762cda8db20f5316449e1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Float_ &gt;</type>
      <name>residuals</name>
      <anchorfile>structscran__variances_1_1FitVarianceTrendResults.html</anchorfile>
      <anchor>a2f01d307daf4f11ad0558e25172e5e79</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_variances::FitVarianceTrendWorkspace</name>
    <filename>structscran__variances_1_1FitVarianceTrendWorkspace.html</filename>
    <templarg>typename Float_</templarg>
  </compound>
  <compound kind="struct">
    <name>scran_variances::ModelGeneVariancesBlockedBuffers</name>
    <filename>structscran__variances_1_1ModelGeneVariancesBlockedBuffers.html</filename>
    <templarg>typename Stat_</templarg>
    <member kind="variable">
      <type>std::vector&lt; ModelGeneVariancesBuffers&lt; Stat_ &gt; &gt;</type>
      <name>per_block</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesBlockedBuffers.html</anchorfile>
      <anchor>a6d216ca738c5f1171c24a1fa5643639c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ModelGeneVariancesBuffers&lt; Stat_ &gt;</type>
      <name>average</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesBlockedBuffers.html</anchorfile>
      <anchor>a2b5059a052625bddab66169fe5b8184e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_variances::ModelGeneVariancesBlockedResults</name>
    <filename>structscran__variances_1_1ModelGeneVariancesBlockedResults.html</filename>
    <templarg>typename Stat_</templarg>
    <member kind="variable">
      <type>std::vector&lt; ModelGeneVariancesResults&lt; Stat_ &gt; &gt;</type>
      <name>per_block</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesBlockedResults.html</anchorfile>
      <anchor>a924b8308d55e630808c783b543ca89c4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>ModelGeneVariancesResults&lt; Stat_ &gt;</type>
      <name>average</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesBlockedResults.html</anchorfile>
      <anchor>ade429098d3a50eb760e69b229d80b49b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_variances::ModelGeneVariancesBuffers</name>
    <filename>structscran__variances_1_1ModelGeneVariancesBuffers.html</filename>
    <templarg>typename Stat_</templarg>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>means</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesBuffers.html</anchorfile>
      <anchor>a2c8bfdbeb8eb86ced05eaeafe7cae484</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>variances</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesBuffers.html</anchorfile>
      <anchor>aab31549e74c37274bccc1a691eecb2a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>fitted</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesBuffers.html</anchorfile>
      <anchor>abbfb61ecb3c2aef9ef4a04efff3bc5c9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Stat_ *</type>
      <name>residuals</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesBuffers.html</anchorfile>
      <anchor>a6db98e3f6ace51821686ceaed178e4d6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_variances::ModelGeneVariancesOptions</name>
    <filename>structscran__variances_1_1ModelGeneVariancesOptions.html</filename>
    <member kind="variable">
      <type>FitVarianceTrendOptions</type>
      <name>fit_variance_trend_options</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesOptions.html</anchorfile>
      <anchor>a2b6b4ecf7ad8a02a930d683382b0e926</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>BlockAveragePolicy</type>
      <name>block_average_policy</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesOptions.html</anchorfile>
      <anchor>aae7a01cbb0c6ffcb0131d11d0c699a7d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>scran_blocks::WeightPolicy</type>
      <name>block_weight_policy</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesOptions.html</anchorfile>
      <anchor>adcc6c3a8d53c6c7b785dedc53f7ef962</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>scran_blocks::VariableWeightParameters</type>
      <name>variable_block_weight_parameters</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesOptions.html</anchorfile>
      <anchor>abe8aea0b4b820b47f52ed7c302d47db0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>block_quantile</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesOptions.html</anchorfile>
      <anchor>a432cca04ce2e5d5b0c897c2c70b5f326</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesOptions.html</anchorfile>
      <anchor>ac7673607b1436c8edde68a3df08519aa</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_variances::ModelGeneVariancesResults</name>
    <filename>structscran__variances_1_1ModelGeneVariancesResults.html</filename>
    <templarg>typename Stat_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>means</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesResults.html</anchorfile>
      <anchor>a36eee93c18e635f5bd3531dd1af129fd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>variances</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesResults.html</anchorfile>
      <anchor>a57c63b5f006d82e5de2ff1a7f4f07871</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>fitted</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesResults.html</anchorfile>
      <anchor>a299a3d92ce9b0936c71209f39b117048</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Stat_ &gt;</type>
      <name>residuals</name>
      <anchorfile>structscran__variances_1_1ModelGeneVariancesResults.html</anchorfile>
      <anchor>a399ae7f3838bb65fed1ce059d3b6572b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>scran_variances</name>
    <filename>namespacescran__variances.html</filename>
    <class kind="struct">scran_variances::ChooseHighlyVariableGenesOptions</class>
    <class kind="struct">scran_variances::FitVarianceTrendOptions</class>
    <class kind="struct">scran_variances::FitVarianceTrendResults</class>
    <class kind="struct">scran_variances::FitVarianceTrendWorkspace</class>
    <class kind="struct">scran_variances::ModelGeneVariancesBlockedBuffers</class>
    <class kind="struct">scran_variances::ModelGeneVariancesBlockedResults</class>
    <class kind="struct">scran_variances::ModelGeneVariancesBuffers</class>
    <class kind="struct">scran_variances::ModelGeneVariancesOptions</class>
    <class kind="struct">scran_variances::ModelGeneVariancesResults</class>
    <member kind="enumeration">
      <type></type>
      <name>BlockAveragePolicy</name>
      <anchorfile>namespacescran__variances.html</anchorfile>
      <anchor>affb6ebee93bfa1a418f68e3265179feb</anchor>
      <arglist></arglist>
      <enumvalue file="namespacescran__variances.html" anchor="affb6ebee93bfa1a418f68e3265179feba4ea6d1161ea24d7599365f574aff6610">MEAN</enumvalue>
      <enumvalue file="namespacescran__variances.html" anchor="affb6ebee93bfa1a418f68e3265179feba4ce64f35e5957a6f3ecc9895dc1bb90f">QUANTILE</enumvalue>
      <enumvalue file="namespacescran__variances.html" anchor="affb6ebee93bfa1a418f68e3265179febab50339a10e1de285ac99d4c3990b8693">NONE</enumvalue>
    </member>
    <member kind="function">
      <type>void</type>
      <name>choose_highly_variable_genes</name>
      <anchorfile>namespacescran__variances.html</anchorfile>
      <anchor>a836ab33e2b96134e9b3ad47df42aeb5e</anchor>
      <arglist>(const std::size_t n, const Stat_ *const statistic, Bool_ *const output, const ChooseHighlyVariableGenesOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Bool_ &gt;</type>
      <name>choose_highly_variable_genes</name>
      <anchorfile>namespacescran__variances.html</anchorfile>
      <anchor>ad2e99f0413753c20e33c217e273ee120</anchor>
      <arglist>(const std::size_t n, const Stat_ *const statistic, const ChooseHighlyVariableGenesOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Index_ &gt;</type>
      <name>choose_highly_variable_genes_index</name>
      <anchorfile>namespacescran__variances.html</anchorfile>
      <anchor>a526febe6eb07121cdd6649a938992ad2</anchor>
      <arglist>(const Index_ n, const Stat_ *const statistic, const ChooseHighlyVariableGenesOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>fit_variance_trend</name>
      <anchorfile>namespacescran__variances.html</anchorfile>
      <anchor>a774f5059db93cedb24f461201f3131fd</anchor>
      <arglist>(const std::size_t n, const Float_ *const mean, const Float_ *const variance, Float_ *const fitted, Float_ *const residuals, FitVarianceTrendWorkspace&lt; Float_ &gt; &amp;workspace, const FitVarianceTrendOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>FitVarianceTrendResults&lt; Float_ &gt;</type>
      <name>fit_variance_trend</name>
      <anchorfile>namespacescran__variances.html</anchorfile>
      <anchor>ada30bbfa19a0c28c23fb69f61355ea54</anchor>
      <arglist>(const std::size_t n, const Float_ *const mean, const Float_ *const variance, const FitVarianceTrendOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>model_gene_variances_blocked</name>
      <anchorfile>namespacescran__variances.html</anchorfile>
      <anchor>a34ffafc531b113e60ae57b3b5f8256be</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;mat, const Block_ *const block, const ModelGeneVariancesBlockedBuffers&lt; Stat_ &gt; &amp;buffers, const ModelGeneVariancesOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>model_gene_variances</name>
      <anchorfile>namespacescran__variances.html</anchorfile>
      <anchor>af48891a8979a6e46b37aaf7c9bc027cb</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;mat, ModelGeneVariancesBuffers&lt; Stat_ &gt; buffers, const ModelGeneVariancesOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>ModelGeneVariancesResults&lt; Stat_ &gt;</type>
      <name>model_gene_variances</name>
      <anchorfile>namespacescran__variances.html</anchorfile>
      <anchor>ab7571d9531d24fb47fa0874ff5039725</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;mat, const ModelGeneVariancesOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>ModelGeneVariancesBlockedResults&lt; Stat_ &gt;</type>
      <name>model_gene_variances_blocked</name>
      <anchorfile>namespacescran__variances.html</anchorfile>
      <anchor>af13889165f096e5af5db48a6046e2d27</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; &amp;mat, const Block_ *const block, const ModelGeneVariancesOptions &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Model per-gene variance in expression</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
