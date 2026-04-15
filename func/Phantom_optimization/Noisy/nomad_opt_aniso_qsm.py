import PyNomad as nomad
import matlab.engine
import os
import sys
import json
import time

# --------------------------- Globals ---------------------------
eng = None
iter_folder = None
txt_file_path = None
counter = 0
best_obj_value = float("inf")
best_params = None


def configure_matlab():
    global eng
    eng = matlab.engine.start_matlab()

    # --- Add your MATLAB toolboxes / code ---
    # FANSI + your ANISO_QSM + dependencies (dipole_kernel_fansi, gradient_calc, load_nii, save_nii, etc.)
    eng.addpath(eng.genpath("Path to FANSI-toolbox-master"), nargout=0)
    eng.addpath(eng.genpath("path to ANISO folder"), nargout=0)
    eng.addpath(eng.genpath("Path to qMRLab"), nargout=0)
    # Folder where you saved:
    #   aniso_nomad_eval.m
    #   aniso_nomad_save_best.m
    eng.addpath("/Path to nomad folder", nargout=0)

    print("MATLAB engine started and paths added.")


def configure_experiment_run(test_fn):
    global iter_folder, txt_file_path

    out_root = "Output Path"
    os.makedirs(out_root, exist_ok=True)

    iter_folder = os.path.join(out_root, test_fn)
    if os.path.exists(iter_folder) and len(os.listdir(iter_folder)) > 0:
        print("Folder already exists and is not empty. Delete it or choose a different name.")
        sys.exit(1)
    os.makedirs(iter_folder, exist_ok=True)

    txt_file_path = os.path.join(out_root, f"{test_fn}.txt")
    with open(txt_file_path, "w") as f:
        f.write("Optimization results (ANISO-QSM NOMAD)\n")

    print("Experiment folder created:", iter_folder)


def log_best_solution(obj_value, iteration, mu1, mu2, lam, it_plus, it_minus, rmse_p, rmse_m, tv_ratio):
    global best_obj_value, best_params

    if obj_value <= best_obj_value:
        best_obj_value = obj_value
        best_params = {
            "mu1": float(mu1),
            "mu2": float(mu2),
            "lambda": float(lam),               # lambda = alpha/mu1
            "alpha1": float(lam * mu1),         # recovered alpha for convenience
            "maxOuterIter": int(it_plus),
            "saveMinusIter": int(it_minus),
            "rmse_xp": float(rmse_p),
            "rmse_xm": float(rmse_m),
            "tv_ratio": float(tv_ratio),
            "objective": float(obj_value),
            "best_iteration": int(iteration),
        }

        with open(txt_file_path, "a") as f:
            f.write(
                f"BEST @ iter {iteration}: OBJ={obj_value:.6g} | "
                f"mu1={mu1:.6g}, mu2={mu2:.6g}, lambda={lam:.6g}, alpha={lam*mu1:.6g}, "
                f"maxOuterIter={it_plus}, saveMinusIter={it_minus} | "
                f"RMSE_Xp={rmse_p:.6g}, RMSE_Xm={rmse_m:.6g}, TVratio={tv_ratio:.6g}\n"
            )

        print("New best solution:", best_params)


def aniso_optimizer(x):
    global counter

    mu1 = x.get_coord(0)
    mu2 = x.get_coord(1)
    lam = x.get_coord(2)  # lambda = alpha/mu1

    # FIX iterations (as you requested)
    it_plus = 100
    it_minus = 10

    # Call MATLAB objective (NO SAVING)
    # Returns: obj, rmse_p, rmse_m, tv_ratio
    obj, rmse_p, rmse_m, tv_ratio = eng.aniso_nomad_eval(
        float(mu1), float(mu2), float(lam), int(it_plus), int(it_minus),
        nargout=4
    )

    obj = float(obj)
    rmse_p = float(rmse_p)
    rmse_m = float(rmse_m)
    tv_ratio = float(tv_ratio)

    print(
        f"Iter {counter}: mu1={mu1:.6g}, mu2={mu2:.6g}, lambda={lam:.6g}, alpha={lam*mu1:.6g} | "
        f"maxOuterIter={it_plus}, saveMinusIter={it_minus} | "
        f"RMSE_Xp={rmse_p:.6g}, RMSE_Xm={rmse_m:.6g}, TVratio={tv_ratio:.6g} | OBJ={obj:.6g}"
    )

    # Log per-iteration numeric sidecar (small, no maps)
    sidecar = {
        "iteration": int(counter),
        "mu1": float(mu1),
        "mu2": float(mu2),
        "lambda": float(lam),
        "alpha1": float(lam * mu1),
        "maxOuterIter": int(it_plus),
        "saveMinusIter": int(it_minus),
        "rmse_xp": float(rmse_p),
        "rmse_xm": float(rmse_m),
        "tv_ratio": float(tv_ratio),
        "objective": float(obj),
    }
    with open(os.path.join(iter_folder, f"iter_{counter:04d}.json"), "w") as jf:
        json.dump(sidecar, jf, indent=2)

    # Update global best (in memory)
    log_best_solution(obj, counter, mu1, mu2, lam, it_plus, it_minus, rmse_p, rmse_m, tv_ratio)

    # Return objective to NOMAD
    x.setBBO(str(obj).encode("UTF-8"))

    counter += 1
    return 1


if __name__ == "__main__":
    # ------------------- CONFIG -------------------
    configure_matlab()
    configure_experiment_run("ANISO_RMSE_ROImean_NOMAD")

    start_time = time.time()

    # ------------------- NOMAD SETUP -------------------
    # Variables: [mu1 (R), mu2 (R), lambda (R)] where lambda = alpha/mu1
    nomad_params = [
        "DIMENSION 3",
        "BB_INPUT_TYPE (R R R)",
        "BB_OUTPUT_TYPE OBJ",
        "MAX_BB_EVAL 150",
        "DISPLAY_DEGREE 1",
        "DISPLAY_ALL_EVAL false",
        "DISPLAY_STATS BBE OBJ",
        "VNS_MADS_SEARCH true",
        "VNS_MADS_SEARCH_TRIGGER 0.75"
    ]

    # Initial guess 
    x0 = [1, 0.7, 0.001]  # mu1, mu2, lambda

    # Bounds (edit if you want)
    lb = [1, 0.1, 1e-6]
    ub = [1e6,  1, 1e-1]

    # ------------------- RUN NOMAD -------------------
    result = nomad.optimize(aniso_optimizer, x0, lb, ub, nomad_params)

    # ------------------- PRINT NOMAD OUTPUT -------------------
    fmt = ["{} = {}".format(n, v) for (n, v) in result.items()]
    print("\nNOMAD results\n" + "\n".join(fmt) + "\n")

    # ------------------- SAVE BEST MAPS ONCE -------------------
    if best_params is None:
        print("No best params found (unexpected).")
        sys.exit(1)

    best_out_dir = os.path.join(iter_folder, "BEST_MAPS")
    os.makedirs(best_out_dir, exist_ok=True)

    # FIX iterations for the final run too
    it_plus = 100
    it_minus = 10

    eng.aniso_nomad_save_best(
        float(best_params["mu1"]),
        float(best_params["mu2"]),
        float(best_params["lambda"]),
        int(it_plus),
        int(it_minus),
        best_out_dir,
        nargout=0
    )

    with open(os.path.join(best_out_dir, "best_params.json"), "w") as jf:
        json.dump(best_params, jf, indent=2)

    end_time = time.time()
    print(f"\nSaved BEST maps in: {best_out_dir}")
    print(f"Best params: {best_params}")
    print(f"Optimization complete in: {end_time - start_time:.3f} seconds")