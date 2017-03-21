// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME Kazu_AdrienDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
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
#include "HourierTracker.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLHourierTracker_Dictionary();
   static void larlitecLcLHourierTracker_TClassManip(TClass*);
   static void *new_larlitecLcLHourierTracker(void *p = 0);
   static void *newArray_larlitecLcLHourierTracker(Long_t size, void *p);
   static void delete_larlitecLcLHourierTracker(void *p);
   static void deleteArray_larlitecLcLHourierTracker(void *p);
   static void destruct_larlitecLcLHourierTracker(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::HourierTracker*)
   {
      ::larlite::HourierTracker *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::HourierTracker));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::HourierTracker", "HourierTracker.h", 30,
                  typeid(::larlite::HourierTracker), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &larlitecLcLHourierTracker_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::HourierTracker) );
      instance.SetNew(&new_larlitecLcLHourierTracker);
      instance.SetNewArray(&newArray_larlitecLcLHourierTracker);
      instance.SetDelete(&delete_larlitecLcLHourierTracker);
      instance.SetDeleteArray(&deleteArray_larlitecLcLHourierTracker);
      instance.SetDestructor(&destruct_larlitecLcLHourierTracker);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::HourierTracker*)
   {
      return GenerateInitInstanceLocal((::larlite::HourierTracker*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::HourierTracker*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLHourierTracker_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::HourierTracker*)0x0)->GetClass();
      larlitecLcLHourierTracker_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLHourierTracker_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLHourierTracker(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::HourierTracker : new ::larlite::HourierTracker;
   }
   static void *newArray_larlitecLcLHourierTracker(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::larlite::HourierTracker[nElements] : new ::larlite::HourierTracker[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLHourierTracker(void *p) {
      delete ((::larlite::HourierTracker*)p);
   }
   static void deleteArray_larlitecLcLHourierTracker(void *p) {
      delete [] ((::larlite::HourierTracker*)p);
   }
   static void destruct_larlitecLcLHourierTracker(void *p) {
      typedef ::larlite::HourierTracker current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::HourierTracker

namespace {
  void TriggerDictionaryInitialization_libKazu_Adrien_Impl() {
    static const char* headers[] = {
"HourierTracker.h",
0
    };
    static const char* includePaths[] = {
"/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/mylarlite/core",
"/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/myLArCV/build/include",
"/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7",
"/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7",
"/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/core/include",
"/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/build/include",
"/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/build/include/../app",
"/Applications/root_v6.08.00/include",
"/Applications/root_v6.08.00/include",
"/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/mylarlite/UserDev/Chimera_Adrien/Adrien/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libKazu_Adrien dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$HourierTracker.h")))  HourierTracker;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libKazu_Adrien dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "HourierTracker.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::HourierTracker", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libKazu_Adrien",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libKazu_Adrien_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libKazu_Adrien_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libKazu_Adrien() {
  TriggerDictionaryInitialization_libKazu_Adrien_Impl();
}
