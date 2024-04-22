add_executable(busyring EXCLUDE_FROM_ALL ring.cpp pm_counters.cpp)
add_dependencies(examples busyring)

find_package(fmt REQUIRED)

target_link_libraries(busyring PRIVATE arbor arborio arborenv arbor-sup ${json_library_name} fmt::fmt-header-only)
