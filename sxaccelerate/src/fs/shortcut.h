// Dummy shortcut.h header for Windows build
#ifndef SHORTCUT_H
#define SHORTCUT_H

#ifdef __cplusplus
extern "C" {
#endif

// Stub functions for Windows shortcut handling
// These functions are not implemented in this build configuration
inline int sxnumlibs_isShortcut(const char* path) {
    // Return 0 (false) - no shortcut support in this build
    (void)path;
    return 0;
}

inline int sxnumlibs_getShortcutTarget(const char* lnkSrc, char** pLnkDest) {
    // Return error code - not implemented
    (void)lnkSrc;
    (void)pLnkDest;
    return -1;
}

#ifdef __cplusplus
}
#endif

#endif // SHORTCUT_H
