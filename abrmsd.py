from ab_rmsd import calc_ab_rmsd
import os
import pandas as pd
import argparse
from DockQ.DockQ import calc_DockQ

def parse_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", type=str, default='rmsd',help='calculate the rmsd or dockQ.')
    parser.add_argument('--native_dir', type=str, default='/user/taosheng/pzz/antibody_data/benchmark_data/igfold_benchmark/pair/')
    parser.add_argument('--pred_dir', type=str, default='/user/taosheng/pzz/antibody_data/predict/igfold_benchmark/alphafold2-m/pair')
    parser.add_argument('--out_path', type=str, default='./igfold.csv')
    args = parser.parse_args()
    return args



def run(
    calc_fn,
    native_dir,
    pred_dir,
    out_path = './result.csv',
    
):
    res = {}
    errors = []
    for pdb_file in os.listdir(pred_dir):
        try:
            id = pdb_file.split('.')[0]
            native_path = os.path.join(native_dir, pdb_file)
            pred_path = os.path.join(pred_dir, pdb_file)
            rmsd = calc_fn(pred_path,native_path)
            res[id] = rmsd
        except Exception as e:
            error_msg = f"[Error] {e}. \n native path: {native_path} \n pred path: {pred_path}"
            print(error_msg)
            errors.append(error_msg)
        
            
    # save error massage
    out_dir = os.path.dirname(out_path)
    with open(os.path.join(out_dir, 'errors.log'), 'w') as f:
        f.write('\n'.join(errors))
        
    # save rmsd statistics
    df = pd.DataFrame(res).T
    df.loc['mean'] = df.mean()
    df.to_csv(out_path,float_format='%.3f')
    pd.options.display.float_format = '{:.2f}'.format
    print(df)
    
def dockQ_batch(
    native_dir,
    pred_dir,
    out_path = './result.csv',
):
    run(calc_DockQ,native_dir,pred_dir,out_path)
    
def rmsd_batch(
    native_dir,
    pred_dir,
    out_path = './result.csv',
):
    run(calc_ab_rmsd,native_dir,pred_dir,out_path)
    
if __name__ == '__main__':
    
    args = parse_arg()
    if args.mode == 'rmsd':
        rmsd_batch(args.native_dir, args.pred_dir, args.out_path)
    else:
        dockQ_batch(args.native_dir, args.pred_dir, args.out_path)
   