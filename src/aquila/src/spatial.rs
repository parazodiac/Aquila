use std::ffi::CStr;
use std::fs::File;
use std::os::raw::c_char;

use std::io::BufWriter;
use std::path::Path;

#[no_mangle]
pub extern "C" fn morans_i(
    weight_file_char: *const c_char,
    values_file_char: *const c_char,
    out_path: *const c_char,
    num_threads: *const c_char,
) -> bool {
    let weight_file_str = unsafe {
        assert!(!weight_file_char.is_null());
        CStr::from_ptr(weight_file_char)
    };

    let values_file_str = unsafe {
        assert!(!values_file_char.is_null());
        CStr::from_ptr(values_file_char)
    };

    let out_file_str = unsafe {
        assert!(!out_path.is_null());
        CStr::from_ptr(out_path)
    };
    let out_file =
        BufWriter::new(File::create(out_file_str.to_str().unwrap()).expect("can't open file"));

    let num_threads_str = unsafe {
        assert!(!num_threads.is_null());
        CStr::from_ptr(num_threads)
    };
    let num_threads: usize = num_threads_str.to_str().unwrap().parse::<usize>().unwrap();

    println!("Using {} threads", num_threads);
    indus::spatial::generate_stats(
        Path::new(weight_file_str.to_str().unwrap()).to_path_buf(),
        Path::new(values_file_str.to_str().unwrap()).to_path_buf(),
        out_file,
        Some("Moransi"),
        num_threads as usize,
    )
    .unwrap();

    true
}

#[no_mangle]
pub extern "C" fn gearys_c(
    weight_file_char: *const c_char,
    values_file_char: *const c_char,
    out_path: *const c_char,
    num_threads: *const c_char,
) -> bool {
    let weight_file_str = unsafe {
        assert!(!weight_file_char.is_null());
        CStr::from_ptr(weight_file_char)
    };

    let values_file_str = unsafe {
        assert!(!values_file_char.is_null());
        CStr::from_ptr(values_file_char)
    };

    let out_file_str = unsafe {
        assert!(!out_path.is_null());
        CStr::from_ptr(out_path)
    };
    let out_file =
        BufWriter::new(File::create(out_file_str.to_str().unwrap()).expect("can't open file"));

    let num_threads_str = unsafe {
        assert!(!num_threads.is_null());
        CStr::from_ptr(num_threads)
    };
    let num_threads: usize = num_threads_str.to_str().unwrap().parse::<usize>().unwrap();

    println!("Using {} threads", num_threads);
    indus::spatial::generate_stats(
        Path::new(weight_file_str.to_str().unwrap()).to_path_buf(),
        Path::new(values_file_str.to_str().unwrap()).to_path_buf(),
        out_file,
        Some("Gearyc"),
        num_threads as usize,
    )
    .unwrap();

    true
}
