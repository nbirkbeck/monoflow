cc_library(
    name = "autorecon-flow",
    srcs = [
        "conversions.cc",
        "filters.cc",
        "flow.cc",
        "flow_variational.cc",
        "multigrid.cc",
        "simplesequence.cc",
        "utils.cc",
    ],
    hdrs = [
        "conversions.h",
        "depthmesh.h",
        "filters.h",
        "flow.h",
        "flow_variational.h",
        "multigrid.h",
        "simplesequence.h",
        "utils.h",
    ],
    copts = ["-fpermissive -DGL_GLEXT_PROTOTYPES"],
    visibility = ["//visibility:public"],
    deps = [
        "//autorecon",
    ],
)
