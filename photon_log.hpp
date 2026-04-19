#pragma once
/**
 * photon_log_cxx11.hpp — Header-only "Drehbuch"-Logger (C++11 kompatibel)
 *
 * - Zwei Overloads für logEnd(): mit/ohne Pose
 *
 * CSV/JSONL bleiben identisch:
 * photon_id, seq, event, x, y, z, dx, dy, dz, end_reason, note
 */

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

struct Pose {
  double x;
  double y;
  double z;
  double dx;
  double dy;
  double dz;
  Pose()
      : x(std::numeric_limits<double>::quiet_NaN()),
        y(std::numeric_limits<double>::quiet_NaN()),
        z(std::numeric_limits<double>::quiet_NaN()),
        dx(std::numeric_limits<double>::quiet_NaN()),
        dy(std::numeric_limits<double>::quiet_NaN()),
        dz(std::numeric_limits<double>::quiet_NaN()) {}
  Pose(double X, double Y, double Z, double DX, double DY, double DZ)
      : x(X), y(Y), z(Z), dx(DX), dy(DY), dz(DZ) {}
};

enum EventKind {
  InitPhotonPose = 0,
  RefractedIntoMaterial,
  PhotonLeavesSpecimen,
  SubSurfaceScattered,
  FresnelReflPointIN,
  FresnelReflPointENDIN,
  AbsorbedInMaterial,
  ScatteredInMaterial,
  ReflectedOnMaterial,
  ReflectedDirect,
  ENDPHOTON
};

enum EndReason { EndNone = 0, Hit, Gone, Absorbed, MaxBounces, OutOfBounds };

inline const char *to_string(EventKind k) {
  switch (k) {
  case InitPhotonPose:
    return "InitPhotonPose";
  case RefractedIntoMaterial:
    return "RefractedIntoMaterial";
  case PhotonLeavesSpecimen:
    return "PhotonLeavesSpecimen";
  case SubSurfaceScattered:
    return "SubSurfaceScattered";
  case FresnelReflPointIN:
    return "FresnelReflPointIN";
  case FresnelReflPointENDIN:
    return "FresnelReflPointENDIN";
  case AbsorbedInMaterial:
    return "AbsorbedInMaterial";
  case ScatteredInMaterial:
    return "ScatteredInMaterial";
  case ReflectedOnMaterial:
    return "ReflectedOnMaterial";
  case ReflectedDirect:
    return "ReflectedDirect";
  case ENDPHOTON:
    return "ENDPHOTON";
  }
  return "Unknown";
}

inline const char *to_string(EndReason r) {
  switch (r) {
  case EndNone:
    return "";
  case Hit:
    return "Hit";
  case Gone:
    return "Gone";
  case Absorbed:
    return "Absorbed";
  case MaxBounces:
    return "MaxBounces";
  case OutOfBounds:
    return "OutOfBounds";
  }
  return "";
}

struct PhotonEvent {
  EventKind kind;
  bool has_pose;
  Pose pose;            // Valid only if has_pose==true
  EndReason end_reason; // Only useful for ENDPHOTON or hard ends
  std::string note;     // optional
  PhotonEvent()
      : kind(InitPhotonPose), has_pose(false), pose(), end_reason(EndNone),
        note() {}
};

class PhotonScriptLogger {
public:
  PhotonScriptLogger() : enabled_(true) {}

  void setEnabled(bool enabled) {
    std::lock_guard<std::mutex> lk(m_);
    enabled_ = enabled;
  }
  bool isEnabled() const { return enabled_; }

  void beginPhoton(uint64_t photon_id) {
    if (!enabled_)
      return;
    std::lock_guard<std::mutex> lk(m_);
    // creates an empty sequence if one does not already exist
    scripts_[photon_id];
  }

  // Event with pose
  void log(uint64_t photon_id, EventKind kind, const Pose &pose,
           const std::string &note = std::string()) {
    if (!enabled_)
      return;
    std::lock_guard<std::mutex> lk(m_);
    PhotonEvent e;
    e.kind = kind;
    e.has_pose = true;
    e.pose = pose;
    e.end_reason = EndNone;
    e.note = note;
    scripts_[photon_id].push_back(e);
  }

  // Event withou pose
  void log(uint64_t photon_id, EventKind kind,
           const std::string &note = std::string()) {
    if (!enabled_)
      return;
    std::lock_guard<std::mutex> lk(m_);
    PhotonEvent e;
    e.kind = kind;
    e.has_pose = false;
    e.end_reason = EndNone;
    e.note = note;
    scripts_[photon_id].push_back(e);
  }

  // END without pose
  void logEnd(uint64_t photon_id, EndReason reason,
              const std::string &note = std::string()) {
    if (!enabled_)
      return;
    std::lock_guard<std::mutex> lk(m_);
    PhotonEvent e;
    e.kind = ENDPHOTON;
    e.has_pose = false;
    e.end_reason = reason;
    e.note = note;
    scripts_[photon_id].push_back(e);
  }

  // END with pose (optional)
  void logEnd(uint64_t photon_id, EndReason reason, const Pose &pose,
              const std::string &note) {
    if (!enabled_)
      return;
    std::lock_guard<std::mutex> lk(m_);
    PhotonEvent e;
    e.kind = ENDPHOTON;
    e.has_pose = true;
    e.pose = pose;
    e.end_reason = reason;
    e.note = note;
    scripts_[photon_id].push_back(e);
  }

  bool writeCSV(const std::string &path) const {
    if (!enabled_)
      return false;
    std::ofstream out(path.c_str());
    if (!out)
      return false;
    out << "photon_id,seq,event,x,y,z,dx,dy,dz,end_reason,note\n";
    std::lock_guard<std::mutex> lk(m_);
    for (typename std::map<uint64_t, std::vector<PhotonEvent>>::const_iterator
             it = scripts_.begin();
         it != scripts_.end(); ++it) {
      uint64_t pid = it->first;
      const std::vector<PhotonEvent> &events = it->second;
      for (size_t i = 0; i < events.size(); ++i) {
        const PhotonEvent &e = events[i];
        out << pid << ',' << i << ',' << to_string(e.kind) << ',';
        if (e.has_pose) {
          out << to_csv_num(e.pose.x) << ',' << to_csv_num(e.pose.y) << ','
              << to_csv_num(e.pose.z) << ',' << to_csv_num(e.pose.dx) << ','
              << to_csv_num(e.pose.dy) << ',' << to_csv_num(e.pose.dz) << ',';
        } else {
          out << ",,,,,,";
        }
        out << to_string(e.end_reason) << ',' << escape_csv(e.note) << '\n';
      }
    }
    return true;
  }

  bool writeJSONL(const std::string &path) const {
    if (!enabled_)
      return false;
    std::ofstream out(path.c_str());
    if (!out)
      return false;
    std::lock_guard<std::mutex> lk(m_);
    for (typename std::map<uint64_t, std::vector<PhotonEvent>>::const_iterator
             it = scripts_.begin();
         it != scripts_.end(); ++it) {
      uint64_t pid = it->first;
      const std::vector<PhotonEvent> &events = it->second;
      for (size_t i = 0; i < events.size(); ++i) {
        const PhotonEvent &e = events[i];
        out << "{";
        out << "\"photon_id\":" << pid << ",";
        out << "\"seq\":" << i << ",";
        out << "\"event\":\"" << to_string(e.kind) << "\",";
        if (e.has_pose) {
          out << "\"pose\":{"
              << "\"x\":" << num_or_null(e.pose.x) << ","
              << "\"y\":" << num_or_null(e.pose.y) << ","
              << "\"z\":" << num_or_null(e.pose.z) << ","
              << "\"dx\":" << num_or_null(e.pose.dx) << ","
              << "\"dy\":" << num_or_null(e.pose.dy) << ","
              << "\"dz\":" << num_or_null(e.pose.dz) << "},";
        } else {
          out << "\"pose\":null,";
        }
        const char *ereas = to_string(e.end_reason);
        if (ereas[0] == '\0')
          out << "\"end_reason\":null,";
        else
          out << "\"end_reason\":\"" << ereas << "\",";
        out << "\"note\":" << json_string(e.note);
        out << "}\n";
      }
    }
    return true;
  }

  std::vector<PhotonEvent> getScript(uint64_t photon_id) const {
    std::lock_guard<std::mutex> lk(m_);
    std::map<uint64_t, std::vector<PhotonEvent>>::const_iterator it =
        scripts_.find(photon_id);
    if (it == scripts_.end())
      return std::vector<PhotonEvent>();
    return it->second;
  }

  void clear() {
    std::lock_guard<std::mutex> lk(m_);
    scripts_.clear();
  }

private:
  static std::string to_csv_num(double v) {
    if (std::isfinite(v)) {
      std::ostringstream oss;
      oss << std::setprecision(17) << v;
      return oss.str();
    }
    return "";
  }

  static std::string escape_csv(const std::string &s) {
    if (s.find_first_of(",\"\n") == std::string::npos)
      return s;
    std::string out = "\"";
    for (size_t i = 0; i < s.size(); ++i) {
      char c = s[i];
      if (c == '"')
        out += "\"\"";
      else
        out += c;
    }
    out += "\"";
    return out;
  }

  static std::string json_string(const std::string &s) {
    std::ostringstream oss;
    oss << "\"";
    for (size_t i = 0; i < s.size(); ++i) {
      char c = s[i];
      switch (c) {
      case '\"':
        oss << "\\\"";
        break;
      case '\\':
        oss << "\\\\";
        break;
      case '\b':
        oss << "\\b";
        break;
      case '\f':
        oss << "\\f";
        break;
      case '\n':
        oss << "\\n";
        break;
      case '\r':
        oss << "\\r";
        break;
      case '\t':
        oss << "\\t";
        break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          // Simple masking for control characters
          oss << "\\u00";
          static const char hex[] = "0123456789ABCDEF";
          unsigned char uc = static_cast<unsigned char>(c);
          oss << hex[(uc >> 4) & 0xF] << hex[uc & 0xF];
        } else {
          oss << c;
        }
      }
    }
    oss << "\"";
    return oss.str();
  }

  static std::string num_or_null(double v) {
    if (std::isfinite(v)) {
      std::ostringstream oss;
      oss << std::setprecision(17) << v;
      return oss.str();
    }
    return "null";
  }

private:
  mutable std::mutex m_;
  bool enabled_;
  std::map<uint64_t, std::vector<PhotonEvent>> scripts_;
};
