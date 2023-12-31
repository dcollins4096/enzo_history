/*
 * Copyright (c) 1999 Apple Computer, Inc. All rights reserved.
 *
 * @APPLE_LICENSE_HEADER_START@
 * 
 * This file contains Original Code and/or Modifications of Original Code
 * as defined in and that are subject to the Apple Public Source License
 * Version 2.0 (the 'License'). You may not use this file except in
 * compliance with the License. Please obtain a copy of the License at
 * http://www.opensource.apple.com/apsl/ and read it before using this
 * file.
 * 
 * The Original Code and all software distributed under the License are
 * distributed on an 'AS IS' basis, WITHOUT WARRANTY OF ANY KIND, EITHER
 * EXPRESS OR IMPLIED, AND APPLE HEREBY DISCLAIMS ALL SUCH WARRANTIES,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, QUIET ENJOYMENT OR NON-INFRINGEMENT.
 * Please see the License for the specific language governing rights and
 * limitations under the License.
 * 
 * @APPLE_LICENSE_HEADER_END@
 */
#ifdef OSX10_4
#include "execinfo_local.h"

/*	Bertrand from vmutils -> CF -> System */

#import <libc.h>
#import <pthread.h>
#import <mach/mach.h>
#include <mach/vm_statistics.h>
#import <malloc/malloc.h>
#import <stdlib.h>

extern void spin_lock(int *);
extern void spin_unlock(int *);
extern void thread_stack_pcs(vm_address_t *, unsigned, unsigned *);

static inline void *allocate_pages(unsigned) __attribute__((always_inline));
static inline void *allocate_pages(unsigned bytes) {
    void *address;
    if (vm_allocate(mach_task_self(), (vm_address_t *)&address, bytes, 
                    VM_MAKE_TAG(VM_MEMORY_ANALYSIS_TOOL)| TRUE)) {
	malloc_printf("*** out of memory while stack logging\n");
	abort();
    } 
    return (void *)address;
}

static inline void deallocate_pages(void *, unsigned) __attribute__((always_inline));
static inline void deallocate_pages(void *ptr, unsigned bytes) {
    vm_deallocate(mach_task_self(), (vm_address_t)ptr, bytes);
}

static inline void copy_pages(const void *, void *, unsigned) __attribute__((always_inline));
static inline void copy_pages(const void *source, void *dest, unsigned bytes) {
    if (vm_copy(mach_task_self(), (vm_address_t)source, bytes, (vm_address_t)dest)) memmove(dest, source, bytes);
}

/***************	Uniquing stack		***********/

#define MAX_COLLIDE	8

#define MAX_NUM_PC	512

static int enter_pair_in_table(unsigned *table, unsigned numPages, unsigned *uniquedParent, unsigned thisPC) {
    // uniquedParent is in-out; return 1 is collisions max not exceeded
    unsigned	base = numPages * vm_page_size / (sizeof(int)*2*2);
    unsigned	hash = base + (((*uniquedParent) << 4) ^ (thisPC >> 2)) % (base - 1); // modulo odd number for hashing
    unsigned	collisions = MAX_COLLIDE;
    while (collisions--) {
        unsigned	*head = table + hash*2;
        if (! head[0] && !head[1]) {
            /* end of chain; store this entry! */
            /* Note that we need to test for both head[0] and head[1] as (0, -1) is a valid entry */
            head[0] = thisPC;
            head[1] = *uniquedParent;
            *uniquedParent = hash;
            return 1;
        }
        if ((head[0] == thisPC) && (head[1] == *uniquedParent)) {
            /* we found the proper entry, the value for the pair is the entry offset */
            *uniquedParent = hash;
            return 1;
        }
        hash++;
        if (hash == base*2) hash = base;
    }
    return 0;
}

unsigned stack_logging_get_unique_stack(unsigned **table, unsigned *table_num_pages, unsigned *stack_entries, unsigned count, unsigned num_hot_to_skip) {
    unsigned	uniquedParent = (unsigned)-1;
    // we skip the warmest entries that are an artefact of the code
    while (num_hot_to_skip--) {
        if (count > 0) { stack_entries++; count--; }
    }
    while (count--) {
        unsigned	thisPC = stack_entries[count];
        while (!enter_pair_in_table(*table, *table_num_pages, &uniquedParent, thisPC)) {
            unsigned	*newTable;
            unsigned	oldBytes = (*table_num_pages) * vm_page_size;
            newTable = allocate_pages(oldBytes*2);
            copy_pages(*table, newTable, oldBytes);
            deallocate_pages(*table, oldBytes);
            *table_num_pages *= 2;
            *table = newTable;
        }
    }
    return uniquedParent;
}

/***************	Logging stack and arguments		***********/

stack_logging_record_list_t *stack_logging_the_record_list = NULL;

int stack_logging_enable_logging = 0;

int stack_logging_dontcompact = 0;

static int stack_logging_spin_lock = 0;

static stack_logging_record_list_t *GrowLogRecords(stack_logging_record_list_t *records, unsigned desiredNumRecords) {
    stack_logging_record_list_t	*new_records;
    unsigned	old_size = records->overall_num_bytes;
    if (desiredNumRecords*sizeof(stack_logging_record_t)+sizeof(stack_logging_record_list_t) < records->overall_num_bytes) return records;
    records->overall_num_bytes += records->overall_num_bytes + vm_page_size; // in order to always get an even number of pages
    new_records = allocate_pages(records->overall_num_bytes);
    copy_pages(records, new_records, old_size);
    deallocate_pages(records, old_size);
    return new_records;
}

static void prepare_to_log_stack(void) {
    if (!stack_logging_the_record_list) {
        unsigned	totalSize = 4 * vm_page_size;
        stack_logging_the_record_list = allocate_pages(totalSize);
        memset(stack_logging_the_record_list, 0, sizeof(stack_logging_record_list_t));
        stack_logging_the_record_list->overall_num_bytes = totalSize;
        stack_logging_the_record_list->uniquing_table_num_pages = 128;
        stack_logging_the_record_list->uniquing_table = allocate_pages(stack_logging_the_record_list->uniquing_table_num_pages * vm_page_size);
    }
}

void stack_logging_log_stack(unsigned type, unsigned arg1, unsigned arg2, unsigned arg3, unsigned result, unsigned num_hot_to_skip) {
    stack_logging_record_t	*rec;
    if (!stack_logging_enable_logging) return;
    // printf("stack_logging_log_stack 0x%x 0x%x 0x%x 0x%x -> 0x%x\n", type, arg1, arg2, arg3, result);
    if (type & stack_logging_flag_zone) {
        // just process it now and be done with it!
        arg1 = arg2; arg2 = arg3; arg3 = 0; type &= ~stack_logging_flag_zone;
    }
    if (type & stack_logging_flag_calloc) {
        // just process it now and be done with it!
        arg1 *= arg2; arg2 = arg3; arg3 = 0; type &= ~stack_logging_flag_calloc;
    }
    if (type & stack_logging_flag_object) {
        unsigned	*class = (unsigned *)arg1;
        arg1 = arg2 + class[5]; // corresponds to the instance_size field
        arg2 = 0; arg3 = 0; type = stack_logging_type_alloc;
    }
    if (type & stack_logging_flag_cleared) {
        type &= ~stack_logging_flag_cleared;
    }
    if (type & stack_logging_flag_handle) {
        if (stack_logging_type_alloc) {
            if (!result) return;
            stack_logging_log_stack(stack_logging_type_alloc, 0, 0, 0, result, num_hot_to_skip+1);
            stack_logging_log_stack(stack_logging_type_alloc, arg1, 0, 0, *((int *)result), num_hot_to_skip+1);
            return;
        }
        if (stack_logging_type_dealloc) {
            if (!arg1) return;
            stack_logging_log_stack(stack_logging_type_dealloc, *((int *)arg1), 0, 0, 0, num_hot_to_skip+1);
            stack_logging_log_stack(stack_logging_type_dealloc, arg1, 0, 0, 0, num_hot_to_skip+1);
            return;
        }
        fprintf(stderr, "*** Unknown logging type: 0x%x\n", type);
    }
    if (type == stack_logging_flag_set_handle_size) {
        if (!arg1) return;
        // Thanks to a horrible hack, arg3 contains the prvious handle value
        if (arg3 == *((int *)arg1)) return;
        stack_logging_log_stack(stack_logging_type_dealloc, arg3, 0, 0, 0, num_hot_to_skip+1);
        stack_logging_log_stack(stack_logging_type_alloc, arg2, 0, 0, *((int *)arg1), num_hot_to_skip+1);
        return;
    }            
    if (type == (stack_logging_type_dealloc|stack_logging_type_alloc)) {
        if (arg1 == result) return; // realloc had no effect, skipping
        if (!arg1) {
            // realloc(NULL, size) same as malloc(size)
            type = stack_logging_type_alloc; arg1 = arg2; arg2 = arg3; arg3 = 0;
	} else {
            // realloc(arg1, arg2) -> result is same as free(arg1); malloc(arg2) -> result
            stack_logging_log_stack(stack_logging_type_dealloc, arg1, 0, 0, 0, num_hot_to_skip+1);
            stack_logging_log_stack(stack_logging_type_alloc, arg2, 0, 0, result, num_hot_to_skip+1);
            return;
        }
    }
    if (type == stack_logging_type_dealloc) {
        // simple free
        if (!arg1) return; // free(nil)
    }
    prepare_to_log_stack();
    spin_lock(&stack_logging_spin_lock);
    stack_logging_the_record_list = GrowLogRecords(stack_logging_the_record_list, stack_logging_the_record_list->num_records + 1);
    rec = stack_logging_the_record_list->records + stack_logging_the_record_list->num_records;
    // We take care of the common case of alloc-dealloc
    if (!stack_logging_dontcompact && stack_logging_the_record_list->num_records && (type == stack_logging_type_dealloc) && arg1 && ((rec-1)->type == stack_logging_type_alloc) && (arg1 == STACK_LOGGING_DISGUISE((rec-1)->address))) {
        stack_logging_the_record_list->num_records--;
        // printf("Erased previous record in alloc-dealloc sequence\n");
    } else {
        unsigned	stack_entries[MAX_NUM_PC];
        unsigned	count = 0;
        rec->type = type;
        if (type == stack_logging_type_dealloc) {
            rec->argument = 0;
            rec->address = STACK_LOGGING_DISGUISE(arg1); // we disguise the address
        } else if (type == stack_logging_type_alloc) {
            rec->argument = arg1;
            rec->address = STACK_LOGGING_DISGUISE(result); // we disguise the address
        } else {
            rec->argument = arg2;
            rec->address = STACK_LOGGING_DISGUISE(arg1); // we disguise the address
        }
	// printf("Before getting samples  0x%x 0x%x 0x%x 0x%x -> 0x%x\n", type, arg1, arg2, arg3, result);
        thread_stack_pcs(stack_entries, MAX_NUM_PC - 1, &count);
        // We put at the bottom of the stack a marker that denotes the thread (+1 for good measure...)
        stack_entries[count++] = (int)pthread_self() + 1;
        /* now let's unique the sample */    
        // printf("Uniquing 0x%x 0x%x 0x%x 0x%x -> 0x%x\n", type, arg1, arg2, arg3, result);
        rec->uniqued_stack = stack_logging_get_unique_stack(&stack_logging_the_record_list->uniquing_table, &stack_logging_the_record_list->uniquing_table_num_pages, stack_entries, count, num_hot_to_skip+2); // we additionally skip the warmest 2 entries that are an artefact of the code
        stack_logging_the_record_list->num_records++;
    }
    spin_unlock(&stack_logging_spin_lock);
}

static kern_return_t default_reader(task_t task, vm_address_t address, vm_size_t size, void **ptr) {
    *ptr = (void *)address;
    return 0;
}

static kern_return_t get_remote_records(task_t task, memory_reader_t reader, stack_logging_record_list_t **records) {
    // sets records
    vm_address_t	*remote_records_address_ref;
    kern_return_t	err;
    *records = NULL;
    err = reader(task, (vm_address_t)&stack_logging_the_record_list, sizeof(vm_address_t), (void **)&remote_records_address_ref);
    if (err) return err;
    if (!*remote_records_address_ref) {
        // printf("stack_logging: no stack record\n");
        return 0;
    }
    // printf("stack_logging: stack records at %p\n", (void *)(*remote_records_address_ref));
    // printf("stack_logging: reading %d bytes\n", sizeof(stack_logging_record_list_t));
    err = reader(task, *remote_records_address_ref, sizeof(stack_logging_record_list_t), (void **)records); // get the list head
    if (err) return err;
    // printf("stack_logging: overall num bytes = %d\n", records->overall_num_bytes);
    return reader(task, *remote_records_address_ref, (*records)->overall_num_bytes, (void **)records);
}

kern_return_t stack_logging_get_frames(task_t task, memory_reader_t reader, vm_address_t address, vm_address_t *stack_frames_buffer, unsigned max_stack_frames, unsigned *num_frames) {
    stack_logging_record_list_t	*records;
    kern_return_t	err;
    unsigned		index;
    unsigned		disguised = STACK_LOGGING_DISGUISE(address);
    if (!reader) reader = default_reader;
    *num_frames = 0;
    err = get_remote_records(task, reader, &records);
    if (err || !records) return err;
    // printf("stack_logging: %d records\n", records->num_records);
    index = 0;
    while (index < records->num_records) {
        stack_logging_record_t	*record = records->records + index;
        if (record->address == disguised) {
	    return stack_logging_frames_for_uniqued_stack(task, reader, record->uniqued_stack, stack_frames_buffer, max_stack_frames, num_frames);
        }
	index++;
    }
    fprintf(stderr, "*** stack_logging: no record found for 0x%x\n", address);
    return 0;
}

kern_return_t stack_logging_enumerate_records(task_t task, memory_reader_t reader, vm_address_t address, void enumerator(stack_logging_record_t, void *), void *context) {
    stack_logging_record_list_t	*records;
    kern_return_t	err;
    unsigned		index;
    unsigned		disguised = STACK_LOGGING_DISGUISE(address);
    if (!reader) reader = default_reader;
    err = get_remote_records(task, reader, &records);
    if (err || !records) return err;
    // printf("stack_logging: %d records\n", records->num_records);
    index = 0;
    while (index < records->num_records) {
        stack_logging_record_t	*record = records->records + index;
        if (!address || (record->address == disguised)) enumerator(*record, context);
	index++;
    }
    return 0;
}

kern_return_t stack_logging_frames_for_uniqued_stack(task_t task, memory_reader_t reader, unsigned uniqued_stack, vm_address_t *stack_frames_buffer, unsigned max_stack_frames, unsigned *num_frames) {
    stack_logging_record_list_t	*records;
    unsigned		*uniquing_table;
    kern_return_t	err;
    if (!reader) reader = default_reader;
    *num_frames = 0;
    err = get_remote_records(task, reader, &records);
    if (err || !records) return err;
    err = reader(task, (vm_address_t)records->uniquing_table, records->uniquing_table_num_pages * vm_page_size, (void **)&uniquing_table);
    if (err) return err;
    while (max_stack_frames && (uniqued_stack != -1)) {
	unsigned	thisPC;
	if ((uniqued_stack * 2 + 1) * sizeof(unsigned) >= records->uniquing_table_num_pages * vm_page_size) {
	    fprintf(stderr, "*** stack_logging: Invalid uniqued stack 0x%x", uniqued_stack);
	    break;
	}
	thisPC = uniquing_table[uniqued_stack * 2];
	uniqued_stack = uniquing_table[uniqued_stack * 2 + 1];
	if (!thisPC && !uniqued_stack) {
	    // Invalid entry
	    fprintf(stderr, "*** stack_logging: Invalid entry 0x%x", thisPC);
	    break;
	}
	stack_frames_buffer[0] = thisPC;
	stack_frames_buffer++;
	(*num_frames)++;
	max_stack_frames--;
    }
    return 0;
}




/*
 * Copyright (c) 2007 Apple Inc. All rights reserved.
 *
 * @APPLE_LICENSE_HEADER_START@
 * 
 * This file contains Original Code and/or Modifications of Original Code
 * as defined in and that are subject to the Apple Public Source License
 * Version 2.0 (the 'License'). You may not use this file except in
 * compliance with the License. Please obtain a copy of the License at
 * http://www.opensource.apple.com/apsl/ and read it before using this
 * file.
 * 
 * The Original Code and all software distributed under the License are
 * distributed on an 'AS IS' basis, WITHOUT WARRANTY OF ANY KIND, EITHER
 * EXPRESS OR IMPLIED, AND APPLE HEREBY DISCLAIMS ALL SUCH WARRANTIES,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, QUIET ENJOYMENT OR NON-INFRINGEMENT.
 * Please see the License for the specific language governing rights and
 * limitations under the License.
 * 
 * @APPLE_LICENSE_HEADER_END@
 */

#include <mach/vm_types.h>
#include <sys/uio.h>

#include <dlfcn.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int backtrace(void** buffer, int size) {
	extern void _thread_stack_pcs(vm_address_t *buffer, unsigned max, unsigned *nb, unsigned skip);
	unsigned int num_frames;
	_thread_stack_pcs((vm_address_t*)buffer, size, &num_frames, 1);
	while (num_frames >= 1 && buffer[num_frames-1] == NULL) num_frames -= 1;
	return num_frames;
}

#if __LP64__
#define _BACKTRACE_FORMAT "%-4d%-35s 0x%016x %s + %u"
#define _BACKTRACE_FORMAT_SIZE 82
#else
#define _BACKTRACE_FORMAT "%-4d%-35s 0x%08x %s + %u"
#define _BACKTRACE_FORMAT_SIZE 65
#endif


static int _backtrace_snprintf(char* buf, size_t size, int frame, const void* addr, const Dl_info* info) {
	char symbuf[19];
	const char* image = "???";
	const char* symbol = symbuf;

	if (info->dli_fname) {
		image = strrchr(info->dli_fname, '/') + 1;
		if (image == NULL) image = info->dli_fname;
	}
	
	if (info->dli_sname) {
		symbol = info->dli_sname;
	} else {
		snprintf(symbuf, sizeof(symbuf), "0x%x", info->dli_saddr);
	}

	return snprintf(buf, size,
			_BACKTRACE_FORMAT,
			frame,
			image,
			addr,
			symbol,
			addr - info->dli_saddr) + 1;
}

char** backtrace_symbols(void* const* buffer, int size) {
	int i;
	size_t total_bytes;
	char** result;
	char** ptrs;
	intptr_t strs;
	Dl_info* info = calloc(size, sizeof (Dl_info));
	
	if (info == NULL) return NULL;
	
	// Compute the total size for the block that is returned.
	// The block will contain size number of pointers to the
	// symbol descriptions.

	total_bytes = sizeof(char*) * size;
	
	// Plus each symbol description
	for (i = 0 ; i < size; ++i) {
		dladdr(buffer[i], &info[i]);
		total_bytes += _BACKTRACE_FORMAT_SIZE + 1;
		if (info[i].dli_sname) total_bytes += strlen(info[i].dli_sname);
	}
	
	result = (char**)malloc(total_bytes);
	if (result == NULL) {
		free(info);
		return NULL;
	}
	
	// Fill in the array of pointers and append the strings for
	// each symbol description.
	
	ptrs = result;
	strs = ((intptr_t)result) + sizeof(char*) * size;
	
	for (i = 0; i < size; ++i) {
		ptrs[i] = (char*)strs;
		strs += _backtrace_snprintf((char*)strs, total_bytes, i, buffer[i], &info[i]);
	}
	
	free(info);
	
	return result;
}

void backtrace_symbols_fd(void* const* buffer, int size, int fd) {
	int i;
	char buf[BUFSIZ];
	Dl_info info;
	struct iovec iov[2];

	iov[0].iov_base = buf;

	iov[1].iov_base = "\n";
	iov[1].iov_len = 1;

	for (i = 0; i < size; ++i) {
		memset(&info, 0, sizeof(info));
		dladdr(buffer[i], &info);

		iov[0].iov_len = _backtrace_snprintf(buf, sizeof(buf), i, buffer[i], &info);
		
		writev(fd, iov, 2);
	}
}



/*
 * Copyright (c) 1999, 2007 Apple Inc. All rights reserved.
 *
 * @APPLE_LICENSE_HEADER_START@
 * 
 * This file contains Original Code and/or Modifications of Original Code
 * as defined in and that are subject to the Apple Public Source License
 * Version 2.0 (the 'License'). You may not use this file except in
 * compliance with the License. Please obtain a copy of the License at
 * http://www.opensource.apple.com/apsl/ and read it before using this
 * file.
 * 
 * The Original Code and all software distributed under the License are
 * distributed on an 'AS IS' basis, WITHOUT WARRANTY OF ANY KIND, EITHER
 * EXPRESS OR IMPLIED, AND APPLE HEREBY DISCLAIMS ALL SUCH WARRANTIES,
 * INCLUDING WITHOUT LIMITATION, ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, QUIET ENJOYMENT OR NON-INFRINGEMENT.
 * Please see the License for the specific language governing rights and
 * limitations under the License.
 * 
 * @APPLE_LICENSE_HEADER_END@
 */

/*	Bertrand from vmutils -> CF -> System */

#include <pthread.h>
#include <mach/mach.h>
#include <mach/vm_statistics.h>
#include <stdlib.h>

#if defined(__i386__) || defined(__x86_64__)
#define FP_LINK_OFFSET 1
#elif defined(__ppc__) || defined(__ppc64__)
#define FP_LINK_OFFSET 2
#else
#error  ********** Unimplemented architecture
#endif

#define	INSTACK(a)	((a) >= stackbot && (a) <= stacktop)
#if defined(__ppc__) || defined(__ppc64__) || defined(__x86_64__)
#define	ISALIGNED(a)	((((uintptr_t)(a)) & 0xf) == 0)
#elif defined(__i386__)
#define	ISALIGNED(a)	((((uintptr_t)(a)) & 0xf) == 8)
#endif

__private_extern__  __attribute__((noinline))
void
_thread_stack_pcs(vm_address_t *buffer, unsigned max, unsigned *nb, unsigned skip)
{
    void *frame, *next;
    pthread_t self = pthread_self();
    void *stacktop = pthread_get_stackaddr_np(self);
    void *stackbot = stacktop - pthread_get_stacksize_np(self);

    *nb = 0;

    /* make sure return address is never out of bounds */
    stacktop -= (FP_LINK_OFFSET + 1) * sizeof(void *);

    /*
     * The original implementation called the first_frame_address() function,
     * which returned the stack frame pointer.  The problem was that in ppc,
     * it was a leaf function, so no new stack frame was set up with
     * optimization turned on (while a new stack frame was set up without
     * optimization).  We now inline the code to get the stack frame pointer,
     * so we are consistent about the stack frame.
     */
#if defined(__i386__) || defined(__x86_64__)
    frame = __builtin_frame_address(0);
#elif defined(__ppc__) || defined(__ppc64__)
    /* __builtin_frame_address IS BROKEN IN BEAKER: RADAR #2340421 */
    __asm__ volatile("mr %0, r1" : "=r" (frame));
#endif
    if(!INSTACK(frame) || !ISALIGNED(frame))
	return;
#if defined(__ppc__) || defined(__ppc64__)
    /* back up the stack pointer up over the current stack frame */
    next = *(void **)frame;
    if(!INSTACK(next) || !ISALIGNED(next) || next <= frame)
	return;
    frame = next;
#endif
    while (skip--) {
	next = *(void **)frame;
	if(!INSTACK(next) || !ISALIGNED(next) || next <= frame)
	    return;
	frame = next;
    }
    while (max--) {
        buffer[*nb] = *(vm_address_t *)(((void **)frame) + FP_LINK_OFFSET);
        (*nb)++;
	next = *(void **)frame;
	if(!INSTACK(next) || !ISALIGNED(next) || next <= frame)
	    return;
	frame = next;
    }
}

void
thread_stack_pcs(vm_address_t *buffer, unsigned max, unsigned *nb)
{
    _thread_stack_pcs(buffer, max, nb, 0);

    // The following prevents thread_stack_pcs() from getting tail-call-optimized into _thread_stack_pcs() on 64-bit environments,
    // thus making the "number of hot frames to skip" be more predictable, giving more consistent backtraces.
    // See <rdar://problem/5364825> "stack logging: frames keep getting truncated" for why this is necessary.
    __asm__ volatile("");
}


#endif /* OSX10_4 */
