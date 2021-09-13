#ifndef LIB_GLOBAL_H
#define LIB_GLOBAL_H

#define BIOUTILS_NAMESPACE bioutils

#define BIOUTILS_BEGIN_SUB_NAMESPACE(name) namespace BIOUTILS_NAMESPACE { namespace name {
#define BIOUTILS_END_SUB_NAMESPACE(name) } }

#define BIOUTILS_BEGIN_NAMESPACE(name) namespace name {
#define BIOUTILS_END_NAMESPACE(name) }

#endif // LIB_GLOBAL_H