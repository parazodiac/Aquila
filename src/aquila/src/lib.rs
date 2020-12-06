// Import dependencies
extern crate libc;

// Modules are other .rs source files
mod spatial;

// Export functions called by R
pub use spatial::string_from_rust;
