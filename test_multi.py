# -*- coding: utf-8 -*-
import multiprocessing
import time

def func(msg):
    return multiprocessing.current_process().name + '-' + msg

if __name__ == "__main__":
    pool = multiprocessing.Pool(processes=4) # 创建4个进程
    results = []
    for i in range(20):
        msg = "process %d" %(i)
        results.append(pool.apply_async(func, (msg, )))
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    print ("Sub-process(es) done.")
    #print (results)
    a=[]
    for res in results:
	a.append(res.get())
        #print (res.get())
    print(a)
