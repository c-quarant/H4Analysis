// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIafsdIcerndOchdIworkdIcdIcquarantdItestbeamH4MyRepodIH4AnalysisdImacrosdIAveragePulseShape_C_ACLiC_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/afs/cern.ch/work/c/cquarant/testbeamH4MyRepo/H4Analysis/macros/./AveragePulseShape.C"

// Header files passed via #pragma extra_include

namespace {
  void TriggerDictionaryInitialization_AveragePulseShape_C_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./AveragePulseShape.C",
0
    };
    static const char* includePaths[] = {
"/afs/cern.ch/sw/lcg/app/releases/ROOT/6.05.02/x86_64-slc6-gcc49-opt/root/include",
"/afs/cern.ch/sw/lcg/app/releases/ROOT/6.05.02/x86_64-slc6-gcc49-opt/root/etc",
"/afs/cern.ch/sw/lcg/app/releases/ROOT/6.05.02/x86_64-slc6-gcc49-opt/root/include",
"/afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/bin/../lib/gcc/x86_64-unknown-linux-gnu/4.9.3/../../../../include/c++/4.9.3",
"/afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/bin/../lib/gcc/x86_64-unknown-linux-gnu/4.9.3/../../../../include/c++/4.9.3/x86_64-unknown-linux-gnu",
"/afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/bin/../lib/gcc/x86_64-unknown-linux-gnu/4.9.3/../../../../include/c++/4.9.3/backward",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/interpreter/cling/include",
"/afs/cern.ch/sw/lcg/app/releases/ROOT/6.05.02/x86_64-slc6-gcc49-opt/root/etc/cling",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/sources/root_v6.05.02/root-6.05.02",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/include",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/sources/root_v6.05.02/root-6.05.02/graf3d/g3d/inc",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/sources/root_v6.05.02/root-6.05.02/gui/gui/inc",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/sources/root_v6.05.02/root-6.05.02/io/io/inc",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/sources/root_v6.05.02/root-6.05.02/core/base/../textinput/src",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/core/base/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/core/rint/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/core/thread/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/io/io/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/sources/root_v6.05.02/root-6.05.02/hist/hist/inc",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/math/mathcore/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/net/net/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/tree/tree/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/math/matrix/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/hist/hist/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/graf2d/graf/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/graf2d/gpad/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/graf3d/g3d/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/tree/treeplayer/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/graf2d/freetype/freetype-2.3.12/include",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/graf2d/asimage/libAfterImage",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/graf2d/asimage/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/hist/histpainter/",
"/afs/cern.ch/work/c/cquarant/testbeamH4MyRepo/H4Analysis/macros/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/sources/root_v6.05.02/root-6.05.02/interpreter/llvm/src/include",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/interpreter/llvm/src/include",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/sources/root_v6.05.02/root-6.05.02/interpreter/llvm/src/tools/clang/include",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/interpreter/llvm/src/tools/clang/include",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/sources/root_v6.05.02/root-6.05.02/interpreter/cling/include",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/etc/cling/cint",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/core/metautils/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/math/minuit/",
"/mnt/build/jenkins/workspace/root-release-6.04/BUILDTYPE/Release/COMPILER/gcc49/LABEL/slc6/build/root_v6.05.02-cmake/graf2d/postscript/",
"/afs/cern.ch/sw/lcg/app/releases/ROOT/6.05.02/x86_64-slc6-gcc49-opt/root/include",
"/afs/cern.ch/work/c/cquarant/testbeamH4MyRepo/H4Analysis/macros/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "./AveragePulseShape.C"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"AveragePulseShape", payloadCode, "@",
"MeanTimeMCP", payloadCode, "@",
"MeanTimeShift", payloadCode, "@",
"PulseShapes", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("AveragePulseShape_C_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_AveragePulseShape_C_ACLiC_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_AveragePulseShape_C_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_AveragePulseShape_C_ACLiC_dict() {
  TriggerDictionaryInitialization_AveragePulseShape_C_ACLiC_dict_Impl();
}
