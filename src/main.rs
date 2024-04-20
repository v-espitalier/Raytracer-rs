use image::{Rgba, RgbaImage, ImageResult};
use std::sync::atomic::{AtomicUsize, Ordering, ATOMIC_USIZE_INIT};
use std::thread;
use std::time::Duration;


/////////////////////////////////////////////////////////////////////////////////////
// Initial Atari code for the Raytracer - BASIC Atari 8 bits - (Fits in 10 lignes)
// https://bunsen.itch.io/raytrace-movie-atari-8bit-by-d-scott-williamson
/////////////////////////////////////////////////////////////////////////////////////
// 1 GR.31:SE.0,7,2:SE.1,3,8:SE.2,13,15:DI.DI(16):F.I=0TO15:REA.R:DI(I)=R:N.I:A=2
// 2 F.N=0TO191:F.M=0TO159:POK.77,0:X=0:Y=-A/25:Z=A:I=SGN(M-80.5)
// 3 U=(M-80.5)/(40*1.3):V=(N-80.5)/(80*1.3):W=1/SQR(U*U+V*V+1):U=U*W:V=V*W:G=1
// 4 E=X-I:F=Y-I:P=U*E+V*F-W*Z:D=P*P-E*E-F*F-Z*Z+1:G.5+(D<=0)*3
// 5 T=-P-SQR(D):G.6+(T<=0)*2:D.0,24,6,30,36,12,42,18,9,33,3,27,45,21,39,15
// 6 X=X+T*U:Y=Y+T*V:Z=Z-T*W:E=X-I:F=Y-I:G=Z:P=2*(U*E+V*F-W*G)
// 7 U=U-P*E:V=V-P*F:W=W+P*G:I=-I:G.4
// 8 IF V<0THEN P=(Y+2)/V:S=INT(X-U*P)+INT(Z-W*P):S=S-INT(S/2)*2:V=-V*(S/2+.3)+.2
// 9 C.3-INT(((3*16)*SQR(V)+DI((M-INT(M/4)*4)+(N-INT(N/4)*4)*4)/3)/16):PL.M,191-N
// 10 N.M:N.N:A=10*(A=-1)+(A-.1)*(A>-1):G.2
/////////////////////////////////////////////////////////////////////////////////////


// Direct port of atari BASIC code to Rust
fn raytracing_cpu_atari_like(factor_res: u32, img_file_prefix: &String) {
        
    let n_img_iter: u32 = 110;     // Initial Atari program behavior
    //let n_img_iter: u32 = 10;
    //let factor_res: u32 = 1;   // Initial resolution
    //let factor_res: u32 = 3;

    // line 1
    // GR.31: Set resolution 160x192 (system instruction)
    let cx: u32 = 160;
    let cy: u32 = 192;

    // Default color
    let pixel_black = Rgba([0, 0, 0, 255]);
    // SE.0,7,2:   Color 0 = blue / dark
    let pixel_blue = Rgba([0, 0, 255, 255]);
    // SE.1,3,8:   Color 1 = purple
    let pixel_purple = Rgba([180, 80, 96, 255]);
    // SE.2,13,15: Color 2 = yellow
    let pixel_yellow = Rgba([160, 128, 0, 255]);


    // DI.DI(16): Declare constants array, dimension 16
    // F.I=0TO15:       // Feed it
    // REA.R
    // DI(I)=R
    //N.I
    // (line below) D.0,24,6,30,36,12,42,18,9,33,3,27,45,21,39,15
    // Dithering constants (tramage)
    let di: [f32; 16] = [0., 24., 6., 30., 36., 12., 42., 18., 9., 33., 3., 27., 45., 21., 39., 15.];

    //A=2
    let mut a : f32 = 2.;

    for img_index in 0..n_img_iter {
	    let mut image : RgbaImage = RgbaImage::new(cx * 2 * factor_res, cy * factor_res);
        
        for (im_x_fact, im_y_fact, pixel) in image.enumerate_pixels_mut() {
//        for n in 0..191 {              // line 2;  F.N=0TO19
//            for m in 0..159 {            // F.M=0TO159
            let m_fact: u32 = ((im_x_fact as f32) / 2.) as u32;
            let n_fact: u32 = 192 * factor_res - im_y_fact;

            let mf : f32 = (m_fact as f32) / (factor_res as f32);
            let nf : f32 = (n_fact as f32)/ (factor_res as f32);
            let m  : u32 = mf as u32;
            let n  : u32 = nf as u32;

            // POK.77,0  :  system instruction, skipped

            let mut x : f32 = 0.;           // X=0
            let mut y : f32 = - a / 25.;    // Y=-A/25:
            let mut z : f32 = a;           // Z=A
            let mut i : f32 = -1.;          // I=SGN(M-80.5)
            if mf >= 80.5 {i = 1.;}

            let mut u : f32 = (mf - 80.5) / (40. * 1.3);    // line 3; U=(M-80.5)/(40*1.3)
            let mut v : f32 = (nf - 80.5) / (80. * 1.3);    // V=(N-80.5)/(80*1.3)
            let mut w : f32 = 1. / f32::sqrt(u * u + v * v + 1.);  // W=1/SQR(U*U+V*V+1)
            u = u * w;                                    // U=U*W
            v = v * w;                                    // V=V*W
            let mut g : f32 = 1.;                         // G=1

            let mut p: f32 = 0.;
            loop {
                let mut e : f32 = x - i;                      // line 4; E=X-I
                let mut f : f32 = y - i;                      // F=Y-I
                p = u * e + v * f - w * z;                    // P=U*E+V*F-W*Z
                let d : f32 = p * p - e * e - f * f - z * z + 1.; // D=P*P-E*E-F*F-Z*Z+1

                // goto:  G.5+(D<=0)*3   // if d<=0, goto line 8 (break)
                if d <= 0. {break;}

                let mut t : f32 = -p - f32::sqrt(d);          // line 5; T=-P-SQR(D)
                // goto: G.6+(T<=0)*2    // if t<=0, goto line 8 (break)
                if t <= 0. {break;}
                // Dithering constants defined here initially in basic: D.0,24,6,30,36,12,42,18,9,33,3,27,45,21,39,15

                x = x + t * u;  // line 6; X=X+T*U
                y = y + t * v;  // Y=Y+T*V
                z = z - t * w;  // Z=Z-T*W
                e = x - i;      // E=X-I
                f = y - i;      // F=Y-I
                g = z;          // G=Z
                p = 2. * (u * e + v * f - w * g);   // P=2*(U*E+V*F-W*G)

                u = u - p * e;  // line 7; U=U-P*E
                v = v - p * f;  // V=V-P*F
                w = w + p * g;  // W=W+P*G
                i = -i;         // I=-I

            }   // goto: G.4

            if v < 0. {
                p = (y + 2.) / v;      // line 8; IF V<0THEN P=(Y+2)/V

                //let mut s : f32 = f32::floor(x - u * p) + f32::floor(z - w * p);  // S=INT(X-U*P)+INT(Z-W*P)
                let mut s : i32 = (f32::floor(x - u * p) + f32::floor(z - w * p)) as i32;  // S=INT(X-U*P)+INT(Z-W*P)

                // s = s - f32::floor(s / 2.) * 2.;   // S=S-INT(S/2)*2
                s = s % 2; // s - f32::floor(s / 2.) * 2.;   // S=S-INT(S/2)*2    // Chessboard

                v = -v * ((s as f32) / 2. + 0.3) + 0.2;    // V=-V*(S/2+.3)+.2
            }

            // line 9 C.3-INT(((3*16)*SQR(V)+DI((M-INT(M/4)*4)+(N-INT(N/4)*4)*4)/3)/16)
            //let m_diff : usize = usize::try_from(m - ((m / 4)) * 4).unwrap();
            //let m_mod : usize = usize::try_from(m % 4).unwrap();
            let m_mod : usize = usize::try_from(m_fact % 4).unwrap();
            //let n_diff : usize = usize::try_from(n - ((n / 4)) * 4).unwrap();
            //let n_mod : usize = usize::try_from(n % 4).unwrap();
            let n_mod : usize = usize::try_from(n_fact % 4).unwrap();

            let di_index : usize = m_mod + n_mod * 4;
            // Pixel color
            let c = 3 - ( (((3. * 16.) * f32::sqrt(v) + di[di_index] / 3.) / 16.) as u32);

            // PL.M,191-N - system instruction: PLOT(M, 191 - N)
            if c == 0
            {
			    *pixel = pixel_black;
            }
            if c == 1
            {
			    *pixel = pixel_blue;
            }
            else if c == 2
            {
			    *pixel = pixel_purple;
            }
            else if c == 3
            {
			    *pixel = pixel_yellow;
            }

            //} // for m in 0..159 {            // line 10; N.M
        //} // for n in 0..191 {            // N.N  (next n)
        } // for (n, m, pixel) in image.enumerate_pixels_mut() {

        // A=10*(A=-1)+(A-.1)*(A>-1)
        if a <= -1.
            {a = 10.;}
        else
            {a = a - 0.1;}

        if (img_file_prefix.len() > 1) {
            let img_f_name : String = format!("{}{:0>3}.png", img_file_prefix, img_index);
            let res : ImageResult<()> = image.save(&img_f_name);
            println!("Saved image : {}", &img_f_name);
        }
     
    }   // goto: G.2      // for img_index in 0..n_img_iter {
} // fn raytracing_cpu_atari_like


// Internal function for multithreading version
fn raytracing_cpu_single_thread(ref_img_index_min : &i32, ref_img_index_max : &i32, ref_n_img_iter : &i32, ref_factor_res: &i32, img_file_prefix: &String) 
{
    let img_index_min: i32 = *ref_img_index_min;
    let img_index_max: i32 = *ref_img_index_max;
    let n_img_iter: i32        = *ref_n_img_iter;
    let factor_res: i32    = *ref_factor_res;
    //println!("Thread img index min,max : {}, {}", img_index_min, img_index_max);

    // line 1
    // GR.31: Set resolution 160x192 (system instruction)
    let cx: i32 = 160;
    let cy: i32 = 192;
        
    for img_index in img_index_min..img_index_max {
        //println!("Image:{}", img_index);
        // DI.DI(16): Declare constants array, dimension 16
        // F.I=0TO15:       // Feed it
        // REA.R
        // DI(I)=R
        //N.I

        //A=2
        let mut a : f32 = 2. - (img_index as f32) * 0.1 / (n_img_iter as f32) * 110.;
        if a <= -1. {a = a + 11.;}
        //let mut a : f32 = 1.;

        // One iteration = one generated image
	    let mut image : RgbaImage = RgbaImage::new((cx * 2 * factor_res) as u32, (cy * factor_res) as u32);
        
        for (im_x_fact, im_y_fact, pixel) in image.enumerate_pixels_mut() {
//        for n in 0..191 {              // line 2;  F.N=0TO19
//            for m in 0..159 {            // F.M=0TO159
            let m_fact: i32 = ((im_x_fact as f32) / 2.) as i32;
            let n_fact: i32 = 192 * factor_res - (im_y_fact as i32);

            let mf : f32 = (m_fact as f32) / (factor_res as f32);
            let nf : f32 = (n_fact as f32)/ (factor_res as f32);
            let m  : i32 = mf as i32;
            let n  : i32 = nf as i32;

            // POK.77,0  :  system instruction, skipped

            let mut x : f32 = 0.;           // X=0
            let mut y : f32 = - a / 25.;    // Y=-A/25:
            let mut z : f32 = a;           // Z=A
            let mut i : f32 = -1.;          // I=SGN(M-80.5)
            if mf >= 80.5 {i = 1.;}

            let mut u : f32 = (mf - 80.5) / (40. * 1.3);    // line 3; U=(M-80.5)/(40*1.3)
            let mut v : f32 = (nf - 80.5) / (80. * 1.3);    // V=(N-80.5)/(80*1.3)
            let mut w : f32 = 1. / f32::sqrt(u * u + v * v + 1.);  // W=1/SQR(U*U+V*V+1)
            u = u * w;                                    // U=U*W
            v = v * w;                                    // V=V*W
            let mut g : f32 = 1.;                         // G=1

            let mut p: f32 = 0.;
            loop {
                let mut e : f32 = x - i;                      // line 4; E=X-I
                let mut f : f32 = y - i;                      // F=Y-I
                p = u * e + v * f - w * z;                    // P=U*E+V*F-W*Z
                let d : f32 = p * p - e * e - f * f - z * z + 1.; // D=P*P-E*E-F*F-Z*Z+1

                // goto:  G.5+(D<=0)*3   // if d<=0, goto line 8 (break)
                if d <= 0. {break;}

                let mut t : f32 = -p - f32::sqrt(d);          // line 5; T=-P-SQR(D)
                // goto: G.6+(T<=0)*2    // if t<=0, goto line 8 (break)
                if t <= 0. {break;}

                x = x + t * u;  // line 6; X=X+T*U
                y = y + t * v;  // Y=Y+T*V
                z = z - t * w;  // Z=Z-T*W
                e = x - i;      // E=X-I
                f = y - i;      // F=Y-I
                g = z;          // G=Z
                p = 2. * (u * e + v * f - w * g);   // P=2*(U*E+V*F-W*G)

                u = u - p * e;  // line 7; U=U-P*E
                v = v - p * f;  // V=V-P*F
                w = w + p * g;  // W=W+P*G
                i = -i;         // I=-I

            }   // goto: G.4

            if v < 0. {
                p = (y + 2.) / v;      // line 8; IF V<0THEN P=(Y+2)/V

                //let mut s : f32 = f32::floor(x - u * p) + f32::floor(z - w * p);  // S=INT(X-U*P)+INT(Z-W*P)

                let mut s : i32 = (f32::floor(x - u * p) + f32::floor(z - w * p)) as i32;  // S=INT(X-U*P)+INT(Z-W*P)

                // s = s - f32::floor(s / 2.) * 2.;   // S=S-INT(S/2)*2
                s = s % 2; // s - f32::floor(s / 2.) * 2.;   // S=S-INT(S/2)*2    // Chessboard

                v = -v * ((s as f32) / 2. + 0.3) + 0.2;    // V=-V*(S/2+.3)+.2
            }

            // line 9 C.3-INT(((3*16)*SQR(V)+DI((M-INT(M/4)*4)+(N-INT(N/4)*4)*4)/3)/16)
            //let m_diff : usize = usize::try_from(m - ((m / 4)) * 4).unwrap();
            //let m_mod : usize = usize::try_from(m % 4).unwrap();
            //let m_mod : usize = usize::try_from(m_fact % 4).unwrap();
            //let n_diff : usize = usize::try_from(n - ((n / 4)) * 4).unwrap();
            //let n_mod : usize = usize::try_from(n % 4).unwrap();
            //let n_mod : usize = usize::try_from(n_fact % 4).unwrap();

            //let di_index : usize = m_mod + n_mod * 4;
            // Pixel color
            //let mut c: u8 = (255. * v) as u8;
            let mut c: u8 = (255. * f32::sqrt(v)) as u8;

            // PL.M,191-N - system instruction: PLOT(M, 191 - N)
            c = 255 - c;

            if c < 50 {
                *pixel = Rgba([c, 0, c, 255]);
            }
            else if c < 100 {
                *pixel = Rgba([50, c - 50, 50, 255]);
            } 
            else if c < 150 {
                *pixel = Rgba([50 + (c - 100), 50 + (c - 100), 100, 255]);
            }
            else {
                *pixel = Rgba([100, 100 + (c - 150), 100 + (c - 100), 255]);
            }

            //} // for m in 0..159 {            // line 10; N.M
        //} // for n in 0..191 {            // N.N  (next n)
        } // for (n, m, pixel) in image.enumerate_pixels_mut() {

        // A=10*(A=-1)+(A-.1)*(A>-1)
        /*
        if a <= -1.
            {a = 10.;}
        else
            {a = a - 0.1;}
        */
     
        if (img_file_prefix.len() > 1) {
            let img_f_name : String = format!("{}{:0>3}.png", img_file_prefix, img_index);
            let res : ImageResult<()> = image.save(&img_f_name);
            println!("Saved image : {}", &img_f_name);
        }
    }   // goto: G.2      // for img_index in 0..n_img_iter {

} // fn raytracing_single_thread()


// Atari -> Rust port (CPU) with higher resolution, higher number of imgs/sec, improved graphics, multithreading
pub fn raytracing_cpu_multithreading(n_img_iter: i32, factor_res: i32, n_threads: u32, img_file_prefix: &String)
{

    let n_img_per_thread: u32 = ((n_img_iter as f32) / (n_threads as f32)).ceil() as u32;

    static GLOBAL_THREAD_COUNT: AtomicUsize = ATOMIC_USIZE_INIT;
    
    for thread_index in 0..n_threads
    {
        let img_index_min: i32 = (n_img_per_thread * thread_index) as i32;
        let mut img_index_max: i32 = (n_img_per_thread * thread_index + n_img_per_thread) as i32;
        if img_index_max >= n_img_iter {img_index_max = n_img_iter - 1;}
        let img_file_prefix_clone: String = img_file_prefix.clone();

        GLOBAL_THREAD_COUNT.fetch_add(1, Ordering::SeqCst);
        let handle = std::thread::spawn( move ||
        {
            raytracing_cpu_single_thread(&img_index_min, &img_index_max, &n_img_iter, &factor_res, &img_file_prefix_clone);
            GLOBAL_THREAD_COUNT.fetch_sub(1, Ordering::SeqCst);
            std::thread::sleep(std::time::Duration::from_millis(1));
            //println!("Thread has terminated {} / {}", thread_index, n_threads);
        });
        //handle.join();
    }
    
    //let n_sec:  u64 = 120;
    //let n_msec: u64 = n_sec * 1000;
    //std::thread::sleep(std::time::Duration::from_millis(n_msec));
    while GLOBAL_THREAD_COUNT.load(Ordering::SeqCst) != 0 {
        thread::sleep(Duration::from_millis(1)); 
    }

} // pub fn raytracing_cpu_multithreading(n_img_iter: i32, factor_res: i32, n_threads: u32, img_file_prefix: &String)


fn main() {

    if (true) {
        // Plain translation from Atari BASIC to Rust
        let factor_res: u32 = 1;   // Original resolution 160x192 (without doubling the horizontal pixels)
        let img_file_prefix: String = String::from("generated_imgs/img_atari_");
        raytracing_cpu_atari_like(factor_res, &img_file_prefix);
    }

    if (true) {
        // Translation to Rust + img quality improvements + multithreading
        let n_img_iter: i32 = 300; // Increase the images rate (x3)
        let factor_res: i32 = 5;   // Increase the resolution to HD
        let n_threads: u32 = 16;   // Number of threads (should match the nb of cpu core)
        let img_file_prefix: String = String::from("generated_imgs/img_cpu_");
        raytracing_cpu_multithreading(n_img_iter, factor_res, n_threads, &img_file_prefix);
    }
 
}