use sysinfo::{ SystemExt, System };
use log::{ warn, info };

pub fn get_available_memory(memory: f64) -> f64 {
    let mut system = System::new_all();
    system.refresh_memory();
    let sys_memory = (system.available_memory() * 1000) as f64;
    info!("Available memory is {} KB", sys_memory / 1000.0);
    if sys_memory <= 500_000_000.0 {
        panic!("Available memory should be greater than 0.5 GB!");
    }

    let available_memory =
        (if memory == 0.0 {
            sys_memory
        } else {
            info!("Requested memory by the user is {} GigaByte", memory);
            if sys_memory < memory * 1000_000_000.0 {
                panic!("Requested memory is greater than the available memory!");
            }
            memory * 1000_000_000.0
        }) - 500_000_000.0;
    available_memory
}

pub fn get_cpus(requested_threads: usize, paired_read_input: bool) -> (usize, usize) {
    let available_cpus = num_cpus::get();
    let used_cpus = if requested_threads == 0 || requested_threads > available_cpus {
        available_cpus
    } else {
        requested_threads
    };
    info!(
        "Available CPUs: {}. Requested CPUs {}. Used CPUs: {}.",
        available_cpus,
        requested_threads,
        used_cpus
    );
    if requested_threads > available_cpus {
        warn!(
            "Used CPUs are redueced to the available CPUs as the requested number is not available."
        );
    }
    let mut reader_threads = 0;
    if used_cpus > 1 {
        reader_threads += 1;
    }
    if used_cpus > 4 && paired_read_input {
        reader_threads += 1;
    }
    let processing_threads = used_cpus - reader_threads;
    if used_cpus > 1 {
        info!(
            "Multi-threading mode is enabled. Readers: {}, Processors {}.",
            reader_threads,
            processing_threads
        );
    }

    (reader_threads, processing_threads)
}
