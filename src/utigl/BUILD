cc_library(
    name = "glcommon",
    srcs = [
        "fbo.cc",
        "ffont.cc",
        "glcommon.cc",
    ],
    hdrs = [
        "fbo.h",
        "ffont.h",
        "glcommon.h",
        "glwindow.h",
        "smartptr.h",
    ],
    copts = ["-DGL_GLEXT_PROTOTYPES"],
    visibility = ["//visibility:public"],
    deps = [
        "@freetype2//:lib",
    ],
)
