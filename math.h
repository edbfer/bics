#pragma once

#include "complex.h"
#include "field.h"

namespace math
{
	template <typename T>
	complex<T> simpson2d(field<T> f);
}