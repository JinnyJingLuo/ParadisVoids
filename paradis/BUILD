cc_library(
    name = "paradis",
    hdrs = glob(["*.h"]),
    srcs = glob(["*.cpp"]),
    copts = ["-std=c++14"],
    visibility = ["//visibility:public"],
    deps = [":dislocation_base", "@ezmath//:ezmath"],
)

cc_library(
    name = "dislocation_base",
    hdrs = glob(["dislocation_base/*.h"]),
    srcs = glob(["dislocation_base/*.cpp"]),
    copts = ["-std=c++14"],
    visibility = ["//visibility:public"],
    deps = ["@ezmath//:ezmath"],
)
