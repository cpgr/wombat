#pragma once

#include <algorithm>
#include <string>
#include <vector>

/**
 * Helper shared by the VE table UserObjects (VERelPermTableUO,
 * VECapillaryPressureTableUO) for reading CSV tables by column NAME rather than
 * position, so the column order written by an upstream tool does not matter.
 *
 * Defined inline in a header (not in an anonymous namespace in each .C) because
 * MOOSE's unity build concatenates the .C files into a single translation unit,
 * where two same-named anonymous-namespace functions would collide.
 */
namespace VECsvColumnLookup
{
/// Return the index of the first column whose header matches any of @p candidates,
/// or -1 if none of them are present in @p names.
inline int
findNamedColumn(const std::vector<std::string> & names,
                const std::vector<std::string> & candidates)
{
  for (const auto & cand : candidates)
  {
    const auto it = std::find(names.begin(), names.end(), cand);
    if (it != names.end())
      return static_cast<int>(std::distance(names.begin(), it));
  }
  return -1;
}
}
