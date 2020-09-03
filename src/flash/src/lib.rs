extern crate libc;

use std;
use std::ffi::CString;
use std::os::raw::c_char;

#[no_mangle]
pub extern fn string_from_rust() -> *const c_char {
    let s = CString::new("Hello ピカチュウ !").unwrap();
    let p = s.as_ptr();
    std::mem::forget(s);
    p
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
