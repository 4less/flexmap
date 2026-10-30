//! buildcmp <fasta> <threads>: builds the index with `default_build` and with `default_build_par`
//! (flexalign's parameters) and checks they are identical; prints both times.
use flexmap::build::{default_build, default_build_par};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (path, threads): (&str, usize) = (&a[1], a[2].parse().unwrap());
    let t = std::time::Instant::now();
    let (m1, r1, i1) = default_build::<31, 15, 16, 4, 12, 16, 2>(path, 4096).unwrap();
    let one = t.elapsed().as_secs_f64();
    let t = std::time::Instant::now();
    let (m2, r2, i2) = default_build_par::<31, 15, 16, 4, 12, 16, 2>(path, 4096, threads).unwrap();
    let par = t.elapsed().as_secs_f64();
    let keys_same = bytemuck::cast_slice::<_, u8>(&m1.keys.data) == bytemuck::cast_slice::<_, u8>(&m2.keys.data);
    let values_same = m1.values.data.len() == m2.values.data.len()
        && m1.values.data.iter().zip(&m2.values.data).all(|(a, b)| a.0 == b.0);
    println!("one thread {one:.1} s, {threads} threads {par:.1} s ({:.1}x); keys identical {keys_same}, values identical {values_same}, names identical {}",
        one / par, i1 == i2 && r1 == r2);
}
