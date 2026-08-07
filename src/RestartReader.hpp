/*******************************************************************************
 * This file is part of CMacIonize
 * Copyright (C) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CMacIonize is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CMacIonize is distributed in the hope that it will be useful,
 * but WITOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file RestartReader.hpp
 *
 * @brief Restart file reader.
 *
 * @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
 */
#ifndef RESTARTREADER_HPP
#define RESTARTREADER_HPP

/*! @brief Enable this to produce a detailed info file that describes everything
 *  that was read by the reader. */
//#define RESTARTREADER_INFO

#include "Error.hpp"

#include <cstdlib>
#include <fstream>
#include <map>
#include <string>

/**
 * @brief Restart file reader.
 */
class RestartReader {
private:
  /*! @brief Underlying input file. */
  std::ifstream _file;

  /*! @brief Restart filename, used in diagnostics. */
  const std::string _filename;

#ifdef RESTARTREADER_INFO
  /*! @brief Detailed info file describing everything that was read by the
   *  reader. */
  std::ofstream _info_file;
#endif

  /**
   * @brief Resolve the requested restart filename.
   *
   * Setting CMAC_RESTART_USE_BACKUP to a non-zero value changes a canonical
   * restart.dump request to restart.0.back. This is intentionally an
   * environment-only recovery mechanism so it does not alter serialized
   * parameters or the restart format.
   */
  inline static std::string resolve_filename(const std::string &filename) {
    const char *use_backup = std::getenv("CMAC_RESTART_USE_BACKUP");
    if (use_backup == nullptr || use_backup[0] == '\0' ||
        std::string(use_backup) == "0") {
      return filename;
    }
    const std::string suffix = "restart.dump";
    if (filename.size() >= suffix.size() &&
        filename.compare(filename.size() - suffix.size(), suffix.size(),
                         suffix) == 0) {
      return filename.substr(0, filename.size() - suffix.size()) +
             "restart.0.back";
    }
    return filename;
  }

  /** @brief Current byte offset, or -1 if unavailable. */
  inline std::streamoff get_offset() {
    const std::streampos position = _file.tellg();
    return position == std::streampos(-1) ? -1 :
                                           static_cast< std::streamoff >(position);
  }

  /** @brief Number of unread bytes, or -1 if unavailable. */
  inline std::streamoff get_remaining_bytes() {
    const std::streampos current = _file.tellg();
    if (current == std::streampos(-1)) {
      return -1;
    }
    _file.seekg(0, std::ios::end);
    const std::streampos end = _file.tellg();
    _file.seekg(current);
    if (end == std::streampos(-1)) {
      return -1;
    }
    return static_cast< std::streamoff >(end - current);
  }

  /** @brief Read exactly size bytes or fail with a useful diagnostic. */
  inline void read_bytes(char *buffer, const std::streamsize size,
                         const char *description) {
    const std::streamoff offset = get_offset();
    _file.read(buffer, size);
    if (!_file) {
      cmac_error("Unexpected end/corruption while reading %s (%lld bytes) "
                 "from restart file \"%s\" at byte offset %lld.",
                 description, static_cast< long long >(size),
                 _filename.c_str(), static_cast< long long >(offset));
    }
  }

public:
  /**
   * @brief Constructor.
   *
   * @param filename Name of the restart file.
   */
  inline RestartReader(const std::string filename)
      : _file(), _filename(resolve_filename(filename)) {
    _file.open(_filename, std::ios::binary);
    if (!_file.is_open()) {
      cmac_error("Failed to open restart file \"%s\" for reading.",
                 _filename.c_str());
    }
    if (_filename != filename) {
      cmac_warning("Restart recovery requested: reading backup file \"%s\" "
                   "instead of \"%s\".",
                   _filename.c_str(), filename.c_str());
    }

#ifdef RESTARTREADER_INFO
    _info_file.open("restart_reader_info.txt");
#endif
  }

  /**
   * @brief Check whether unread data remains in the restart file.
   *
   * This is useful for optional data appended to the end of a restart file,
   * while retaining compatibility with restart files written by older
   * versions.
   *
   * @return True if at least one unread byte remains.
   */
  inline bool has_more_data() {
    return _file.peek() != std::ifstream::traits_type::eof();
  }

  /**
   * @brief General read function for basic template data types.
   *
   * @return Value read from the restart file.
   */
  template < typename _datatype_ > _datatype_ read() {
    _datatype_ value{};
    read_bytes(reinterpret_cast< char * >(&value), sizeof(_datatype_),
               "restart value");
#ifdef RESTARTREADER_INFO
    _info_file << sizeof(_datatype_) << "\n";
#endif
    return value;
  }
};

/**
 * @brief RestartReader::read specialization for a bool.
 *
 * @return bool read from the restart file.
 */
template <> inline bool RestartReader::read() {
  uint_least8_t boolean = read< uint_least8_t >();
#ifdef RESTARTREADER_INFO
  _info_file << "bool\n";
#endif
  return boolean > 0;
}

/**
 * @brief RestartReader::read specialization for a string.
 *
 * @return std::string read from the restart file.
 */
template <> inline std::string RestartReader::read() {
  const auto size = read< std::string::size_type >();
  const std::streamoff remaining = get_remaining_bytes();
  if (remaining < 0 || size > static_cast< std::string::size_type >(remaining)) {
    cmac_error("Invalid string length %zu in restart file \"%s\" at byte "
               "offset %lld; only %lld bytes remain. The restart is truncated "
               "or incompatible with this build.",
               static_cast< size_t >(size), _filename.c_str(),
               static_cast< long long >(get_offset()),
               static_cast< long long >(remaining));
  }
  std::string string(size, '\0');
  if (size > 0) {
    read_bytes(&string[0], static_cast< std::streamsize >(size),
               "restart string");
  }
#ifdef RESTARTREADER_INFO
  _info_file << "string\n";
#endif
  return string;
}

/**
 * @brief RestartReader::read specialization for a map of strings.
 *
 * @return std::map<std::string, std::string> read from the restart file.
 */
template <> inline std::map< std::string, std::string > RestartReader::read() {
  const auto size = read< std::map< std::string, std::string >::size_type >();
  const std::streamoff remaining = get_remaining_bytes();
  const size_t minimum_bytes_per_entry =
      2 * sizeof(std::string::size_type);
  if (remaining < 0 ||
      size > static_cast< std::map< std::string, std::string >::size_type >(
                 remaining / static_cast< std::streamoff >(
                                 minimum_bytes_per_entry))) {
    cmac_error("Invalid map length %zu in restart file \"%s\" at byte offset "
               "%lld. The restart is truncated or incompatible with this "
               "build.",
               static_cast< size_t >(size), _filename.c_str(),
               static_cast< long long >(get_offset()));
  }
  std::map< std::string, std::string > map;
  for (std::map< std::string, std::string >::size_type i = 0; i < size; ++i) {
    const std::string key = read< std::string >();
    const std::string value = read< std::string >();
    map[key] = value;
  }
#ifdef RESTARTREADER_INFO
  _info_file << "map\n";
#endif
  return map;
}

#endif // RESTARTREADER_HPP
