include(FetchContent)

FetchContent_Declare(
  tatami
  GIT_REPOSITORY https://github.com/tatami-inc/tatami
  GIT_TAG master # ^3.0.0
)

FetchContent_Declare(
  tatami_stats
  GIT_REPOSITORY https://github.com/tatami-inc/tatami_stats
  GIT_TAG master # ^1.1.0
)

FetchContent_Declare(
  WeightedLowess 
  GIT_REPOSITORY https://github.com/LTLA/CppWeightedLowess
  GIT_TAG master # ^2.0.0
)

FetchContent_Declare(
  scran_core_utils
  GIT_REPOSITORY https://github.com/libscran/core_utils
  GIT_TAG master
)

FetchContent_MakeAvailable(tatami)
FetchContent_MakeAvailable(tatami_stats)
FetchContent_MakeAvailable(WeightedLowess)
FetchContent_MakeAvailable(scran_core_utils)
